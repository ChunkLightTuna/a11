#!/usr/bin/env python3
from __future__ import print_function
import math
import cairo
import sys
from subprocess import Popen, PIPE
import time
import threading
from builtins import str, int, ValueError, len, bytes, range, min

__author__ = 'Bence Cserna'
TIME_LIMIT = 60

# Maximum window size
MAX_X = 800
MAX_Y = 800

# University of New Hampshire brand colors
BRAND_DARK = [37, 55, 70]
BRAND_LIGHT1 = (162, 170, 173)
BRAND_LIGHT2 = (214, 210, 196)
BRAND_ORANGE = [247, 122, 5]


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)


def debug(*objs):
    print("DEBUG: ", *objs, file=sys.stdout)


def filter_input(lines):
    return [line for line in lines if not line.startswith('/') and len(line) > 0]


def parse_meta_information(data):
    circles = []
    circle_count = int(data.pop(0).strip().split(" ")[-1])
    for _ in range(0, circle_count):
        x, y, r = map(float, data.pop(0).strip()[3:].split(" "))
        circles.append((x, y, r))

    lines = []
    line_count = int(data.pop(0).strip().split(" ")[-1])
    for _ in range(0, line_count):
        x0, y0, x1, y1 = map(float, data.pop(0).strip()[3:].split(" "))
        lines.append((x0, y0, x1, y1))

    return circles, lines


def read_problem(data):
    first_line = data[0]

    try:
        circles = []
        lines = []

        if first_line.startswith("/ Meta information:"):
            data.pop(0)
            circles, lines = parse_meta_information(data)

        filtered_input = filter_input(data)
        width = int(filtered_input.pop(0).strip())
        height = int(filtered_input.pop(0).strip())

        map = []
        for x in range(0, height):
            line = filtered_input[x].strip()

            row = []
            for char in line:
                if char in '_':
                    row.append(False)
                elif char in '#':
                    row.append(True)
                else:
                    error("Invalid character: %s" % char)

            if not len(row) == width:
                error("Invalid row length %s expected %s" % (len(row), width))
                error(row)
                error(line)
                sys.exit(-1)
            map.append(row)

        if not len(row) == width:
            error("Invalid number of lines %s expected %s" % (len(map), height))
            error(row)
            error(line)
            sys.exit(-1)

        return circles, lines, map
    except ValueError:
        error("Invalid input!")
        sys.exit(0)


def read_solution(data):
    circles = []
    circle_count = int(data.pop(0).strip().split(" ")[-1])
    for _ in range(0, circle_count):
        x, y, r = map(float, data.pop(0).strip().split(" "))
        circles.append((x, y, r))

    lines = []
    line_count = int(data.pop(0).strip().split(" ")[-1])
    for _ in range(0, line_count):
        x0, y0, x1, y1 = map(float, data.pop(0).strip().split(" "))
        lines.append((x0, y0, x1, y1))

    return circles, lines


def print_usage():
    print("Usage: executable")


def execute_application(executable, params, process, raw_input_line):
    cmd = [" ".join(["./" + executable])]

    debug("Command: " + str(cmd))
    process[0] = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    start = time.time()
    debug("Waiting for termination")
    out, err = process[0].communicate(str.encode(raw_input_line))

    end = time.time()
    debug("Execution time: %s seconds" % (end - start))

    raw_policy = bytes.decode(out).split("\n")

    return read_solution(raw_policy)


def test_application(executable, args, process):
    raw_input = sys.stdin.readlines()
    raw_input_line = ''.join(raw_input)

    original_circles, original_lines, map = read_problem(raw_input)
    circles, lines = execute_application(executable, None, process, raw_input_line)  # Handle termination

    graphics = Graphics(map, (original_circles, original_lines), (circles, lines))
    graphics.draw("ransac-%s" % time.strftime("%Y%m%d-%H%M%S"))


class Validator:
    def __init__(self, executable, args):
        self.args = args
        self.executable = executable
        self.process = [None]

    def run(self, timeout):
        def validate(executable, args, process):
            test_application(executable, args, process)

        thread = threading.Thread(target=validate,
                                  kwargs={"executable": self.executable, "args": self.args, "process": self.process})
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
            self.process[0].terminate()  # This is not thread safe!
            error("Time limit reached: %s seconds" % timeout)
            thread.join()


class Graphics:
    def __init__(self, point_map, original, generated):
        self.generated = generated
        self.original = original
        self.point_map = point_map
        self.width = len(point_map[0])
        self.height = len(point_map)
        self.cell_size = int(min(MAX_X / self.width, MAX_Y / self.height))
        self.offset = self.cell_size / 2

    def draw_points(self, ctx):
        for w in range(0, self.width):
            for h in range(0, self.height):
                x0 = self.cell_size * w
                y0 = self.cell_size * h

                ctx.rectangle(x0, y0, self.cell_size, self.cell_size)

                if self.point_map[h][w]:
                    ctx.fill()
                else:
                    ctx.stroke()

    def draw_circles(self, ctx, circles):
        ctx.set_line_width(0.9)
        for circle in circles:
            self.draw_circle(ctx, circle[0], circle[1], circle[2])

    def draw_circle(self, ctx, x, y, r):
        ctx.move_to(x * self.cell_size + self.offset, y * self.cell_size + self.offset)
        ctx.arc(x * self.cell_size + self.offset, y * self.cell_size + self.offset,
                self.cell_size * r,
                0, 2 * math.pi)
        ctx.stroke()

    def draw_lines(self, ctx, lines):
        for line in lines:
            self.draw_line(ctx, line[0], line[1], line[2], line[3])

    def draw_line(self, ctx, x0, y0, x1, y1):
        ctx.set_line_width(0.9)
        ctx.move_to(x0 * self.cell_size + self.offset, y0 * self.cell_size + self.offset)
        ctx.line_to(x1 * self.cell_size + self.offset, y1 * self.cell_size + self.offset)
        ctx.stroke()

    def draw_grid(self, ctx):
        for w in range(0, self.width):
            for h in range(0, self.height):
                ctx.set_line_width(0.1)
                ctx.set_source_rgba(BRAND_ORANGE[0] / 255, BRAND_ORANGE[1] / 255, BRAND_ORANGE[2] / 255, 0.5)

                x0 = self.cell_size * w
                y0 = self.cell_size * h

                ctx.rectangle(x0, y0, x0 + self.cell_size, y0 + self.cell_size)
                ctx.stroke()

    def mirror_circles_y(self, circles):
        return [(circle[0], self.height - circle[1] - 1, circle[2]) for circle in circles]

    def mirror_lines_y(self, lines):
        return [(line[0], self.height - line[1] - 1, line[2], self.height - line[3] - 1) for line in lines]

    def draw(self, file_name):
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width * self.cell_size, self.height * self.cell_size)
        ctx = cairo.Context(surface)

        # paint background
        ctx.set_source_rgb(BRAND_DARK[0] / 255, BRAND_DARK[1] / 255, BRAND_DARK[2] / 255)
        ctx.rectangle(0, 0, self.width * self.cell_size, self.height * self.cell_size)
        ctx.fill()

        ctx.set_line_width(0.1)
        ctx.set_source_rgba(BRAND_ORANGE[0] / 255, BRAND_ORANGE[1] / 255, BRAND_ORANGE[2] / 255, 0.5)
        self.draw_points(ctx)

        # Draw reference lines and circles
        ctx.set_source_rgba(0, 1, 0, 0.5)  # Green
        self.draw_circles(ctx, self.original[0])
        self.draw_lines(ctx, self.original[1])

        # ctx.set_source_rgb(0, 0, 1)
        # self.draw_circles(ctx, self.generated[0])
        # self.draw_lines(ctx, self.generated[1])

        # Draw given lines and circles
        ctx.set_source_rgb(1, 1, 1)  # White
        self.draw_circles(ctx, self.mirror_circles_y(self.generated[0]))
        self.draw_lines(ctx, self.mirror_lines_y(self.generated[1]))

        surface.write_to_png("%s.png" % file_name)
        debug("Image generated: %s.png" % file_name)


def main():
    if len(sys.argv) < 2:
        print_usage()
        return

    try:
        executable, args = sys.argv[1], sys.argv[2:]
    except ValueError:
        error(["Invalid input: "] + sys.argv)
        print_usage()
        return

    validator = Validator(executable, args)
    validator.run(TIME_LIMIT)


if __name__ == "__main__":
    main()
