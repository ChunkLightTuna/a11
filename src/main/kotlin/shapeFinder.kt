import java.io.BufferedReader
import java.io.FileInputStream
import java.io.InputStream
import java.security.SecureRandom
import java.util.*
import kotlin.system.exitProcess

class shapeFinder {
    companion object {
        //magic numbers
        val SUFFICIENT_NUMBER_OF_INLIERS = 7
        val LARGE_ENOUGH_GAP_BETWEEN_POINTS = 3.5
        val SUFFICIENT_DENSITY_OF_INLIERS = .7
        val RANSAC_ATTEMPTS = 500
        val ERROR = 4

        var DEBUG = false

        @JvmStatic fun main(args: Array<String>) {

            if (args.isNotEmpty() && args[0].endsWith(".w"))
                System.setIn(FileInputStream(args[0]))
            DEBUG = (args.size > 1 && args[1] == "--DEBUG")

            val inputStream = System.`in`

            val shapeFinder = shapeFinder()
            val image = shapeFinder.parse(inputStream)
            inputStream.close()


//            if (DEBUG) {
//                println("${closePoints.size} out of ${image.pixels.size}")
//                image.DEBUG(closePoints)
//            }

            val lines = ransac(image)

            println("number of circles: 0")
            println("number of lines: ${lines.size}")

            lines.forEach { println("$it") }

        }

        private fun ransac(image: Image): Collection<LineSegment> {
            val out = mutableListOf<LineSegment>()

            var lewps = 0
            while (lewps < RANSAC_ATTEMPTS && image.pixels.size > SUFFICIENT_NUMBER_OF_INLIERS) {
                lewps += 1

                var line = image.bootStrapLine()
                var closePoints = image.closePoints(line)

                var i = 0
                var found = closePoints.size
                while (i < RANSAC_ATTEMPTS) {
                    if (DEBUG)
                        println("#$lewps: found ${closePoints.size} points out of ${image.pixels.size}")

                    line = Line.fit(closePoints)
                    closePoints = image.closePoints(line)

                    i += if (found == closePoints.size || found <= 2) RANSAC_ATTEMPTS else 1

                    found = closePoints.size
                }
                if (closePoints.size < SUFFICIENT_NUMBER_OF_INLIERS) {
                    continue
                }

                val pointsSorted = closePoints.sortedBy { it -> line.projection(it).x }

                var last = pointsSorted.first()

                val listOfLists = mutableListOf(mutableListOf(last))

                pointsSorted.forEach {
                    if (it.distance(last) > LARGE_ENOUGH_GAP_BETWEEN_POINTS) {
                        listOfLists.add(mutableListOf())
                    }
                    listOfLists.last().add(it)

                    last = it
                }

                if (DEBUG) {
                    listOfLists.forEach {
                        println("List:$it")
                    }
                }

                val best = listOfLists.filter { list ->
                    list.size > SUFFICIENT_NUMBER_OF_INLIERS &&
                            list.size / (list.first().distance(list.last())) > SUFFICIENT_DENSITY_OF_INLIERS
                }.maxBy { it -> it.size } ?: continue

                val segment = LineSegment(best.first(), best.last())

                best.forEach {
                    image.pixels.remove(it)
                }

                out.add(segment)
            }
            return out
        }
    }

    fun parse(inputStream: InputStream): Image {
        val lines: BufferedReader = inputStream.bufferedReader()

        var line = lines.readLine()
        while (line.startsWith('/')) {
            line = lines.readLine()
        }
        val cols = line.toInt()
        val rows = lines.readLine().toInt()

        val pixels = hashSetOf<Point>()
        val grid = Array(rows, { BooleanArray(cols) })

        (0 until rows).forEach { y ->
            val next = lines.readLine()
            next.forEachIndexed { x, c ->
                if (c == '#') {
                    pixels.add(Point(x, y))
                }
                grid[x][y] = (c == '#')
            }
        }

        return Image(cols, rows, /*grid,*/ pixels)


    }
}

class Image(val cols: Int, val rows: Int, /*val grid: Array<BooleanArray>,*/ val pixels: HashSet<Point>, val random: Random = Random()) {
    fun bootStrapLine(): Line {
        assert(pixels.size > 1)
        val i1 = random.nextInt(pixels.size)
        var i2 = random.nextInt(pixels.size)
        while (i1 == i2)
            i2 = random.nextInt(pixels.size)

        return LineSegment(Pair(pixels.elementAt(i1), pixels.elementAt(i2))).toLine()
    }

    fun closePoints(line: Line): List<Point> {
        return pixels.filter { line.distance(it) < shapeFinder.ERROR }
    }

    fun debug(points: Collection<Point>) {

        val grid = Array(rows, { i -> Array(cols, { i -> false }).toBooleanArray() })

        points.forEach { grid[it.x.toInt()][it.y.toInt()] = true }

        print(' ')
        (0 until cols).forEach { print('_') }
        println()
        grid.forEach {
            print('|')
            print(it.map { if (it) 'o' else ' ' }.toCharArray())
            print("|\n")
        }
        print(' ')
        (0 until cols).forEach { print('_') }
        println()
    }
}

data class Point(val x: Double, val y: Double) {
    constructor(x: Int, y: Int) : this(x.toDouble(), y.toDouble())

    fun distance(other: Point): Double {
        return Math.sqrt(Math.pow((other.x - x), 2.0) + Math.pow((other.y - y), 2.0))
    }
}

data class Line(val m: Double, val b: Double, val denom: Double = Math.sqrt(m * m + 1)) {
    fun distance(p: Point): Double {
        return Math.abs(m * p.x - p.y + b) / denom
    }

    fun projection(point: Point): Point {
        val x = point.x + m / (Math.pow(m, 2.0) + 1) * (point.y - b - m * point.x)
        return Point(x, m * x + b)
    }

    companion object {
        fun fit(points: Collection<Point>): Line {
            val xBar = points.map(Point::x).average()
            val yBar = points.map(Point::y).average()
            val sizeMinusOne = points.size - 1

            val Sxx = points.map { it -> Math.pow(it.x - xBar, 2.0) }.sum() / sizeMinusOne
            val Sxy = points.map { it -> (it.x - xBar) * (it.y - yBar) }.sum() / sizeMinusOne
            val Syy = points.map { it -> Math.pow(it.y - yBar, 2.0) }.sum() / sizeMinusOne

            //Note that βˆ1 is the slope(m) of the line you want and βˆ0 is its y-intercept(b).
            val β1 =
                    (Syy - Sxx + Math.sqrt(Math.pow(Syy - Sxx, 2.0) + 4 * Math.pow(Sxy, 2.0))) /
                            (2 * Sxy)
            val β0 = yBar - β1 * xBar

            val line = Line(β1, β0)


            val starHats = points.map { it ->
                val x = it.x + β1 / (Math.pow(β1, 2.0) + 1) * (it.y - β0 - β1 * it.x)
                Point(x, β1 * x + β0)
            }

            //fuck yea works.
            starHats.forEach { it ->
                assert(line.distance(it) < .00001)
            }

            return line
        }
    }
}

data class LineSegment(val a: Point, val b: Point) {
    constructor(pair: Pair<Point, Point>) : this(pair.first, pair.second)

    fun Double.format(): String? = java.lang.String.format("%.${1}f", this)

    fun toLine(): Line {
        val m = b.y - a.y / b.x - a.x
        val b = a.y - m * a.x
        return Line(m, b)
    }

    override fun toString(): String {
        return "${a.x.format()} ${a.y.format()} ${b.x.format()} ${b.y.format()}"
    }
}
