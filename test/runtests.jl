using GeohashHilbert
using Random
using Test

const GHH = GeohashHilbert

randlat() = -90 + 180 * rand()
randlon() = -180 + 360 * rand()
randprec() = rand(1:20)
randgeohash(n, bits_per_char = 2) = String(rand('0':'3', n))

# converting between int/string and back should reproduce starting values
function test_int_str_conversion()
    for _ in 1:20
        x = rand(0,5000)
        bpc = rand([2,4,6])
        base = 2^bits_per_char
        chars_needed = floor(Int, 1 + log(base, x))
        nchar = rand(chars_needed : (chars_needed + 5))
        x_str = GHH.int_to_str(x, nchar, bpc)
        @test length(x_str) == nchar
        x_prime = GHH.str_to_int(x_str, bpc)
        @test x == x_prime
    end
    return nothing
end

# converting between xy and integer should be lossless
function test_int_xy_conversion()
    for _ in 1:20
        k = rand(2:20)
        n = 2^k
        int = rand(0:(n^2 - 1))
        x, y = GHH.int_to_xy(int, n)
        int_prime = GHH.xy_to_int(x, y, n)
        @test int_prime == int
    end
    return nothing
end

# For lat/lons which are centers of geohash cells, encoding then decoding
# should match input coordinates.
# We can find lat/lons which are centers of cells by encoding and then decoding
# random points.
function test_cell_centers_encode_decode(bits_per_char = 2)
    for _ in 1:50
        lon, lat = randlon(), randlat()
        prec = randprec()
        geohash = GHH.encode(lon, lat, prec, bits_per_char)
        cell_center = GHH.decode(geohash, bits_per_char)
        geohash2 = GHH.encode(cell_center..., prec, bits_per_char)
        @test geohash == geohash2
        cell_center2 = GHH.decode(geohash2, bits_per_char)
        @test cell_center2 == cell_center
    end
    return nothing
end

# Make sure we match Python geohash hilbert on encoding, particularly on edges
# and corners of cells.
function test_match_python_encode(bits_per_char = 2)
    # each item is (lon, lat, prec, python encode at 2 bits_per_char)
    # this covers a variety of nasty edge cases, boundaries between cells, 
    # and points in each quadrant of the world (wrt equator and prime meridian)
    llp_g = [
        [90, 47, 2, "22"],
        [90, 45, 2, "22"],
        [0, 90, 6, "211111"],
        [0, -90, 6, "322222"],
        [-180, 47, 5, "11000"],
        [180, 47, 5, "22333"],
        [180, 90, 7, "2222222"],
        [13, 47, 21, "210032131130033230313"],
        [-13, 47, 21, "123301202203300103020"],
        [13, -47, 21, "321103020023322121020"],
        [-13, -47, 21, "012230313310011212313"]
    ]
    for tup in llp_g
        @test GHH.encode(tup[1:3]...) == tup[4]
    end
    return nothing
end

Random.seed!(47)
test_cell_centers_encode_decode()
test_match_python_encode()
println("Great job!")
