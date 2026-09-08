#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from fractions import Fraction
from pathlib import Path


def parse_resolution(value: str) -> tuple[int, str]:
	match = re.search(r"\d+", value)
	if match is None:
		raise ValueError(f"resolution must contain digits: {value!r}")
	return int(match.group()), value


def convert_count(count: int, old_resolution: int, new_resolution: int) -> int:
	result = Fraction(count - 1, 1) * Fraction(old_resolution, new_resolution) + 1
	if result.denominator != 1:
		raise ValueError(
			f"converted count is not an integer: {count} with {old_resolution} -> {new_resolution}"
		)
	return result.numerator


def update_line(line: str, old_resolution: int, new_resolution: int, old_token: str, new_token: str) -> str:
	updated = line.replace(old_token, new_token)
	updated = updated.replace(f"_a_{old_resolution}_", f"_a_{new_resolution}_")

	nqrs_match = re.match(r"^(\s*btp\((\d+)\)%nqrs\s*=\s*)([^!\n]*?)(\s*,\s*)$", updated)
	if nqrs_match:
		prefix, _, values_text, suffix = nqrs_match.groups()
		values = [part.strip() for part in values_text.split(",") if part.strip()]
		if len(values) != 3:
			raise ValueError(f"expected three nqrs values in line: {line.rstrip()}")
		converted_values = [str(convert_count(int(value), old_resolution, new_resolution)) for value in values]
		return f"{prefix}{',   '.join(converted_values)}{suffix}"

	npml_match = re.match(r"^(\s*btp\((\d+)\)%npml\s*=\s*)(\d+)(\s*[,/]\s*)$", updated)
	if npml_match:
		prefix, _, value_text, suffix = npml_match.groups()
		converted_value = int(Fraction(int(value_text) * old_resolution, new_resolution))
		return f"{prefix}{converted_value}{suffix}"

	return updated


def process_file(path: Path, old_resolution_text: str, new_resolution_text: str) -> Path:
	old_resolution, _ = parse_resolution(old_resolution_text)
	new_resolution, _ = parse_resolution(new_resolution_text)

	original_text = path.read_text()
	updated_lines = [
		update_line(line, old_resolution, new_resolution, old_resolution_text, new_resolution_text)
		for line in original_text.splitlines(keepends=True)
	]
	updated_text = "".join(updated_lines)

	new_path = path.with_name(path.name.replace(old_resolution_text, new_resolution_text))
	if new_path != path:
		path.rename(new_path)
		path = new_path

	if updated_text != original_text:
		path.write_text(updated_text)

	return path


def main(argv: list[str]) -> int:
	if len(argv) != 3:
		print("usage: python change_resolution.py old_resolution new_resolution", file=sys.stderr)
		return 1

	old_resolution_text = argv[1]
	new_resolution_text = argv[2]
	old_resolution, _ = parse_resolution(old_resolution_text)
	new_resolution, _ = parse_resolution(new_resolution_text)

	for path in sorted(Path.cwd().glob("*.in")):
		original_text = path.read_text()
		updated_lines = [
			update_line(line, old_resolution, new_resolution, old_resolution_text, new_resolution_text)
			for line in original_text.splitlines(keepends=True)
		]
		updated_text = "".join(updated_lines)
		new_path = path.with_name(path.name.replace(old_resolution_text, new_resolution_text))

		if new_path != path:
			path.rename(new_path)
			path = new_path

		if updated_text != original_text:
			path.write_text(updated_text)
			print(f"updated {path.name}")
		elif new_path != path:
			print(f"renamed {path.name}")

	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))
