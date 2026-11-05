#!/usr/bin/env python3
"""Extract XYZ structures from a CJKQC JSON Lines database.

Usage:
  JKjsonl2xyz.py DATABASE.jsonl OUTPUT_DIRECTORY

Writes one XYZ file per record that carries a usable geometry and prints
the absolute path of each written file, one per line. Records without a
usable geometry are skipped with a warning. Duplicate structure names are
an error. Exits with status 1 when no structure could be written.
"""

import json
import math
import os
import re
import sys


def fail(message):
  print("error: " + message, file=sys.stderr)
  sys.exit(1)


def safe_name(name, line_number):
  if not isinstance(name, str) or not name.strip():
    return "record_{}".format(line_number)
  name = os.path.basename(name.strip())
  name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
  name = name.strip("._-")
  if not name:
    return "record_{}".format(line_number)
  return name


def record_name(record, line_number):
  name = record.get("name")
  if isinstance(name, str) and name.strip():
    return safe_name(name, line_number)
  source_path = record.get("source_path")
  if isinstance(source_path, str) and source_path.strip():
    stem = os.path.splitext(os.path.basename(source_path.strip()))[0]
    if stem:
      return safe_name(stem, line_number)
  return safe_name(None, line_number)


def extract_geometry(record):
  geometry = record.get("geometry")
  if not isinstance(geometry, dict):
    return None
  symbols = geometry.get("symbols")
  coordinates = geometry.get("coordinates")
  if not isinstance(symbols, list) or not isinstance(coordinates, list):
    return None
  if not symbols or len(coordinates) != 3 * len(symbols):
    return None
  for symbol in symbols:
    if not isinstance(symbol, str) or not symbol:
      return None
  values = []
  for coordinate in coordinates:
    if isinstance(coordinate, bool) or not isinstance(coordinate, (int, float)):
      return None
    coordinate = float(coordinate)
    if not math.isfinite(coordinate):
      return None
    values.append(coordinate)
  return symbols, values


def main(argv):
  if len(argv) != 3 or argv[1] in ("-h", "--help"):
    print("Usage: JKjsonl2xyz.py DATABASE.jsonl OUTPUT_DIRECTORY", file=sys.stderr)
    return 2

  database = argv[1]
  outdir = argv[2]
  if not os.path.isfile(database):
    fail("database does not exist: {}".format(database))
  try:
    os.makedirs(outdir, exist_ok=True)
  except OSError as error:
    fail("cannot create output directory {}: {}".format(outdir, error))

  written = 0
  names = {}
  try:
    handle = open(database, "r", encoding="utf-8")
  except OSError as error:
    fail("cannot read {}: {}".format(database, error))

  with handle:
    for line_number, line in enumerate(handle, 1):
      line = line.strip()
      if not line:
        continue
      try:
        record = json.loads(line)
      except json.JSONDecodeError as error:
        print("warning: {}:{}: cannot parse JSON ({}); skipping".format(
          database, line_number, error), file=sys.stderr)
        continue
      if not isinstance(record, dict):
        print("warning: {}:{}: record is not a JSON object; skipping".format(
          database, line_number), file=sys.stderr)
        continue

      geometry = extract_geometry(record)
      if geometry is None:
        print("warning: {}:{}: no usable geometry; skipping".format(
          database, line_number), file=sys.stderr)
        continue

      name = record_name(record, line_number)
      if name in names:
        fail("duplicate structure name '{}' at lines {} and {}".format(
          name, names[name], line_number))
      names[name] = line_number

      symbols, coordinates = geometry
      path = os.path.abspath(os.path.join(outdir, name + ".xyz"))
      try:
        with open(path, "w", encoding="utf-8") as xyz:
          xyz.write("{}\n".format(len(symbols)))
          xyz.write("{}\n".format(name))
          for atom in range(len(symbols)):
            xyz.write("{} {:.10f} {:.10f} {:.10f}\n".format(
              symbols[atom],
              coordinates[3 * atom],
              coordinates[3 * atom + 1],
              coordinates[3 * atom + 2]))
      except OSError as error:
        fail("cannot write {}: {}".format(path, error))
      print(path)
      written += 1

  if written == 0:
    fail("no usable structures found in {}".format(database))
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv))
