import java.io.BufferedReader
import java.io.FileInputStream
import java.io.InputStream
import java.util.*

fun main(args: Array<String>) {
    if (args.isNotEmpty() && args[0].endsWith(".w"))
        System.setIn(FileInputStream(args[0]))

    val inputStream = System.`in`

    val shapeFinder = shapeFinder()
    val image = shapeFinder.parse(inputStream)
    inputStream.close()

    val lines = shapeFinder.ransac(image)

    println("number of circles: 0")
    println("number of lines: ${lines.size}")

    lines.forEach { println("$it") }
}

fun Double.format(): String? = java.lang.String.format("%.${1}f", this)

class shapeFinder {
    companion object {
        //magic numbers
        val SUFFICIENT_NUMBER_OF_INLIERS = 7
        val LARGE_ENOUGH_GAP_BETWEEN_POINTS = 3.5
        val SUFFICIENT_DENSITY_OF_INLIERS = .7
        val RANSAC_ATTEMPTS = 1000
        val REFIT_ATTEMPTS = 100
    }

    fun ransac(image: Image): Collection<LineSegment> {
        val out = mutableListOf<LineSegment>()

        var lewps = 0
        mainLoop@ while (lewps < RANSAC_ATTEMPTS && image.pixels.size >= SUFFICIENT_NUMBER_OF_INLIERS) {
            lewps += 1

            var fitLine: Line = image.bootStrapLine()
            var matchingPoints: List<Point> = image.closePoints(fitLine)

            var i = 0
            refit@ while (i < REFIT_ATTEMPTS) {

                fitLine = Line.fit(matchingPoints)
                matchingPoints = image.closePoints(fitLine)

                if (matchingPoints.size < 2) {
                    break@refit
                }
                i += 1
            }
            if (matchingPoints.size < SUFFICIENT_NUMBER_OF_INLIERS) {
                continue@mainLoop
            }

//            println("# of refits:$i")

            val pointsSorted = matchingPoints.sortedBy { it -> fitLine.projection(it).x }

            var last = pointsSorted.first()
            val listOfLists = mutableListOf(mutableListOf(last))

            pointsSorted.drop(1).forEach {
                if (it.distance(last) > LARGE_ENOUGH_GAP_BETWEEN_POINTS) {
                    listOfLists.add(mutableListOf())
                }
                listOfLists.last().add(it)

                last = it
            }

            val best = listOfLists.filter { list ->
                list.size >= SUFFICIENT_NUMBER_OF_INLIERS &&
                        list.size / (list.first().distance(list.last())) > SUFFICIENT_DENSITY_OF_INLIERS
            }.maxBy { it -> it.size } ?: continue@mainLoop

            val segment = LineSegment(fitLine.projection(best.first()), fitLine.projection(best.last()))

            best.forEach {
                image.pixels.remove(it)
            }
            out.add(segment)
        }
        return out
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
                    pixels.add(Point(x, rows - y - 1))
                }
                grid[x][y] = (c == '#')
            }
        }

        return Image(pixels)
    }
}

class Image(val pixels: HashSet<Point>, val random: Random = Random()) {
    fun bootStrapLine(): Line {
        assert(pixels.size > 1)
        val i1 = random.nextInt(pixels.size)
        var i2 = random.nextInt(pixels.size)
        while (i1 == i2)
            i2 = random.nextInt(pixels.size)

        return LineSegment(Pair(pixels.elementAt(i1), pixels.elementAt(i2))).toLine()
    }

    fun closePoints(line: Line): List<Point> {
        return pixels.filter { line.distance(it) <= shapeFinder.LARGE_ENOUGH_GAP_BETWEEN_POINTS }
    }
}

class Point(val x: Double, val y: Double) {
    constructor(x: Int, y: Int) : this(x.toDouble(), y.toDouble())

    fun distance(other: Point): Double {
        return Math.sqrt(Math.pow((other.x - x), 2.0) + Math.pow((other.y - y), 2.0))
    }
}

class Line(val m: Double, val b: Double) {

    fun intersect(line: Line): Point {
        val x = (line.b - b) / (m - line.m)
        val y = m * x + b
        return Point(x, y)
    }

    fun distance(p: Point): Double {
        return projection(p).distance(p)
    }

    fun projection(p: Point): Point {
        val x = (p.y + 1 / m * p.x - b) / (m + 1 / m)

        assert(intersect(LineSegment(Point(x, m * x + b), p).toLine()).distance(Point(x, m * x + b)) < .0005)

        return Point(x, m * x + b)
    }

    companion object {
        fun fit(points: Collection<Point>): Line {
            val xBar = points.map(Point::x).average()
            val yBar = points.map(Point::y).average()
            val sizeMinusOne = points.size - 1

            val Sxx = points.fold(0.0, { i, p -> i + Math.pow(p.x - xBar, 2.0) }) / sizeMinusOne
            val Sxy = points.fold(0.0, { i, p -> i + (p.x - xBar) * (p.y - yBar) }) / sizeMinusOne
            val Syy = points.fold(0.0, { i, p -> i + Math.pow(p.y - yBar, 2.0) }) / sizeMinusOne

            //Note that βˆ1 is the slope(m) of the line you want and βˆ0 is its y-intercept(b).
            val β1 = (Syy - Sxx + Math.sqrt(Math.pow(Syy - Sxx, 2.0) + 4 * Math.pow(Sxy, 2.0))) / (2 * Sxy)

            val β0 = yBar - β1 * xBar

            return Line(β1, β0)
        }
    }

    override fun toString(): String {
        return "y = ${m.format()}x ${if (b < 0) '-' else '+'} ${Math.abs(b).format()}"
    }
}

class LineSegment(val a: Point, val b: Point) {
    constructor(pair: Pair<Point, Point>) : this(pair.first, pair.second)

    fun toLine(): Line {
        val m = (b.y - a.y) / (b.x - a.x)
        val b = a.y - (m * a.x)

        return Line(m, b)
    }

    override fun toString(): String {
        return "${a.x.format()} ${a.y.format()} ${b.x.format()} ${b.y.format()}"
    }
}
