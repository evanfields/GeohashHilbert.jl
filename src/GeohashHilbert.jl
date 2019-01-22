module GeohashHilbert

using StaticArrays

"sing Base.ImmutableDict in place of Dict for `ORDERING` and 
`HILBERT_RECURSION` seems to provide a ~2x performance gain. Strangely,
there's not a simple constructor taking a Dict to an ImmutableDict. So,
this function recursively converts `Dict`s to `ImmutableDict`s."
function dict_to_idict(dict::T where T <: Dict)
    ks = keys(dict)
    dict_vt = valtype(dict)
    if dict_vt <: Dict
        new_vt = Base.ImmutableDict{keytype(dict_vt), valtype(dict_vt)}
    else
        new_vt = dict_vt
    end
    idict = Base.ImmutableDict{keytype(dict), new_vt}()
    for k in ks
        idict = Base.ImmutableDict(idict, k => dict_to_idict(dict[k]))
    end
    return idict
end
function dict_to_idict(notdict)
    return notdict
end

# define the ordering of each orientation of a hilbert curve
# note that there are multiple possible global orientations of the Hilbert curve
# corresponding to rotations and reflections
# we use an orientation with global shape Π starting in the lower left
@enum HilbertOrientation UUp ULeft UDown URight # U ] Π [
@enum Quadrant LowerLeft LowerRight UpperLeft UpperRight
const ORDERING = dict_to_idict(Dict(
	UUp => Dict(UpperRight => 1, LowerRight => 2, LowerLeft => 3, UpperLeft => 4),
	ULeft => Dict(LowerLeft => 1, LowerRight => 2, UpperRight => 3, UpperLeft => 4),
	UDown => Dict(LowerLeft => 1, UpperLeft => 2, UpperRight => 3, LowerRight => 4),
	URight => Dict(UpperRight => 1, UpperLeft => 2, LowerLeft => 3, LowerRight => 4)
))
# define the recursive nature of the Hilbert curve
const HILBERT_RECURSION = dict_to_idict(Dict(
	UUp => Dict(UpperRight => URight, LowerRight => UUp, LowerLeft => UUp, UpperLeft => ULeft),
	ULeft => Dict(LowerLeft => UDown, LowerRight => ULeft, UpperRight => ULeft, UpperLeft => UUp),
	UDown => Dict(LowerLeft => ULeft, UpperLeft => UDown, UpperRight => UDown, LowerRight => URight),
	URight => Dict(UpperRight => UUp, UpperLeft => URight, LowerLeft => URight, LowerRight => UDown)
))
# e.g. in a UDown oriented curve, the lower left quadrant is a ULeft oriented curve


# equivalent to converting x to base 4 and then a string
"Convert an integer to a string of specified length and bits per character encoding."
int_to_str(x::Int, nchar, bits_per_char = 2) = string(x; base = 2^bits_per_char, pad = nchar)
str_to_int(s::AbstractString, bits_per_char = 2) = parse(Int, s, base = 2^bits_per_char)

"""
Encode a lon-lat as a geohash. Higher `precision` leads to a finer grained encoding;
in particular the unwrapped lon-lat plane is divided into `4^precision` squares on a side.
Thus the number of bits used is `2 * precision`.
"""
function encode(lon, lat, precision, bits_per_char = 2)
    # only supported for now
    @assert bits_per_char == 2
	@assert -90 <= lat <= 90
	@assert -180 <= lon <= 180
	encode_bits = precision * bits_per_char
	n = 2^(encode_bits ÷ 2)
	x, y = lonlat_to_xy(lon, lat, n)
	curve_spot = xy_to_int(x, y, n)
	return int_to_str(curve_spot, precision)
end

"""
Convert lon-lat coordinates to a (rounded) `x,y` integer point in the
`n` by `n` grid covering lon-lat space. Returned `x,y` are in [1...n].
"""
function lonlat_to_xy(lon, lat, n)
	# awkward 1+floor rounding to keep xy mapping for lon-lat exactly on
	# cell boundaries consistent with python package GeohashHilbert
	x = min(n, 1 + floor(Int, (lon + 180) / 360 * n))
	y = min(n, 1 + floor(Int, (lat + 90) / 180 * n))
	return x, y
end

"""
Convert `x,y` coordinates in the `n` by `n` grid covering lon-lat space
to longitude-latitude coordinates. The coordinates returned correspond to the
center of the `(x,y)`-th grid rectangle.
"""
function xy_to_lonlat(x, y, n)
    lon_grid_size = 360 / n
    lon = -180 + lon_grid_size * (x - .5)
    lat_grid_size = 180 / n
    lat = -90 + lat_grid_size * (y - .5)
    return lon, lat
end

"""
Convert coordinates `(x,y)` to an integer representing location along the Hilbert
curve filling a `n` by `n` grid. `n` must be a power of 2 and `x,y` should be
integers in `[1...n]`.
"""
function xy_to_int(x::Int, y::Int, n::Int)
	@assert 1 <= x <= n
	@assert 1 <= y <= n
	
	cur_orientation = UDown
	cur_quadrant_size = n ÷ 2
    # steps is steps along the curve, starting with the lower left most point
    # being 0 steps along the curve
    # it might be more natural to use 1-based Julian-style indexing, but given
    # the goal of matching the Python package as closely as possible, this 0-based
    # ends up simplifying the overall code.
	steps = 0
	
	while cur_quadrant_size > 0
		quadrant = get_quadrant(x, y, cur_quadrant_size)
		n_previous_quadrants = ORDERING[cur_orientation][quadrant] - 1
		# steps to get to this quadrant
		steps += n_previous_quadrants * cur_quadrant_size^2
		# now iterate within the quadrant
		if x > cur_quadrant_size
			x -= cur_quadrant_size
		end
		if y > cur_quadrant_size
			y -= cur_quadrant_size
		end
		cur_orientation = HILBERT_RECURSION[cur_orientation][quadrant]
		cur_quadrant_size ÷= 2
	end
	
	return steps
end

"""
Convert an integer `t` representing position along the Hilbert curve filling
a `n` by `n` grid (as always, with upside-down U shape, starting in the lower 
left) to `x,y` coordinates in `[1...n]` by `[1...n]` with `(1,1)` representing
the lower left.
"""
function int_to_xy(t::Int, n::Int)
	if !(1 <= t <= n * n)
        error("t passed to int_to_xy must be in [1,n^2]. Got (t,n) = $((t,n))")
    end
	
	cur_orientation = UDown
	cur_quadrant_size = n ÷ 2
	x = 1
	y = 1
	
	while cur_quadrant_size > 0
		pts_per_quadrant = cur_quadrant_size^2
		quadrant_index = 1
        # >= check because t is steps from the first xy square
        # so if eg t is exactly pts_per_quadrant, then the xy to return is
        # the first point of the second quadrant, not the last point of
        # the first quadrant.
		while t >= pts_per_quadrant
			quadrant_index += 1
			t -= pts_per_quadrant
		end
		quadrant = find_quadrant(quadrant_index, cur_orientation)
		if quadrant == UpperLeft || quadrant == UpperRight
			y += cur_quadrant_size
		end
		if quadrant == LowerRight || quadrant == UpperRight
			x += cur_quadrant_size
		end
		cur_orientation = HILBERT_RECURSION[cur_orientation][quadrant]
		cur_quadrant_size ÷= 2
	end
	
	return x, y
end

# for a given quadrant index (1-4) and orientation,
# find the quadrant (eg UpperLeft) of the corresponding index
# there's probably a small absolute but large relative efficiency gain
# to be had by making a REVERSE_ORDERING dict or using a different data structure
function find_quadrant(quad_index, orientation)::Quadrant
	@assert 1 <= quad_index <= 4
	for quad in instances(Quadrant)
		if ORDERING[orientation][quad] == quad_index
			return quad
		end
	end
	error("Couldn't find position of quad index $(quad_index) for orientation $(orientation)")
end
	
@inline function get_quadrant(x, y, quad_size)
	x <= quad_size && y <= quad_size && return LowerLeft
	x > quad_size && y <= quad_size && return LowerRight
	x <= quad_size && y > quad_size && return UpperLeft
	return UpperRight
end

"Given a `geohash` string at a specified `bits_per_char`, return the coordinates of the
corresponding geohash cell's center as a tuple `(lon, lat)`."
function decode(geohash, bits_per_char = 2)
    # currently only 2 bits per char supported
    @assert bits_per_char == 2
    precision = length(geohash)
    curve_spot = str_to_int(geohash)
    n = 2^(precision * bits_per_char ÷ 2)
    x, y = int_to_xy(curve_spot, n)
    return xy_to_lonlat(x, y, n)
end

function decode_exactly(geohash, bits_per_char = 2)
    # currently only 2 bits per char supported
    @assert bits_per_char == 2
    precision = length(geohash)
    lon, lat = decode(geohash, bits_per_char)
    lon_rect_size, lat_rect_size = cell_size_deg(precision, bits_per_char)
    return lon, lat, lon_rect_size / 2, lat_rect_size / 2
end

"Compute the size in longitude/latitude degrees of geohash rectangles
for a given precision level and bits per character (default 2). Return
a tuple `(lon_size, lat_size)`"
function cell_size_deg(precision, bits_per_char = 2)
    n = 2^(precision * bits_per_char ÷ 2)
    return 360 / n, 180 / n
end

function neighbours(geohash, bits_per_char)
    lon, lat, lon_err, lat_err = decode_exactly(geohash, bits_per_char)
    
    north = lat + 2 * lat_err
    south = lat - 2 * lat_err
    east = lon + 2 * lon_err
    west = lon - 2 * lon_err
    
    # wrap around the lon = +/- 180 line
    east > 180 && (east -= 360)
    west < -180 && (west += 360)
    
    neigh_dict = Dict{String, String}(
        "east" => encode(east, lat, prec, bits_per_char),
        "west" => encode(west, lat, prec, bits_per_char)
    )
    
    if north <= 90 # input cell isn't already at the north pole
        neigh_dict = merge(neigh_dict, Dict{String, String}(
            "north-east" => encode(east, north, prec, bits_per_char),
            "north" => encode(lon, north, prec, bits_per_char),
            "north-west" => encode(west, north, prec, bits_per_char)
        ))
    end
    
    if south >= -90 # input cell isn't already at the south pole
        neigh_dict = merge(neigh_dict, Dict{String, String}(
            "south-west" => encode(west, south, prec, bits_per_char),
            "south" => encode(lon, south, prec, bits_per_char),
            "south-east" => encode(east, south, prec, bits_per_char)
        ))
    end
    
    return neigh_dict
    
end

function rectangle(geohash, bits_per_char)

    lon, lat, lon_err, lat_err = decode_exactly(geohash, bits_per_char)
    
    return Dict{Any, Any}(
        "type" => "Feature",
        "properties" => Dict{Any, Any}(
            "code" => geohash,
            "lon" => lon,
            "lat" => lat,
            "lon_err" => lon_err,
            "lat_err" => lat_err,
            "bits_per_char" => bits_per_char,
        ),
        "bbox" => (
            lon - lon_err,  # bottom left
            lat - lat_err,
            lon + lon_err,  # top right
            lat + lat_err,
        ),
        "geometry" => Dict{Any, Any}(
            "type" => "Polygon",
            "coordinates" => [[
                (lon - lon_err, lat - lat_err),
                (lon + lon_err, lat - lat_err),
                (lon + lon_err, lat + lat_err),
                (lon - lon_err, lat + lat_err),
                (lon - lon_err, lat - lat_err),
            ]],
        ),
    )
end

end # module
