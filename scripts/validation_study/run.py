"""Master runner for Peterson's validation studies."""

import argparse
from typing import List

import ex_4_1
import ex_4_3
import ex_4_24
import matplotlib.pyplot as plt

EXAMPLES = {
    "4.1": ex_4_1.run,
    "4.3": ex_4_3.run,
    "4.24": ex_4_24.run,
}


def main():
    parser = argparse.ArgumentParser(description="Run Peterson's validation examples.")
    parser.add_argument(
        "examples",
        nargs="*",
        choices=list(EXAMPLES.keys()) + ["all"],
        default=["all"],
        help="The examples to run (default: all).",
    )

    args = parser.parse_args()

    to_run: List[str] = args.examples
    if "all" in to_run:
        to_run = list(EXAMPLES.keys())

    for ex_name in to_run:
        print(f"\n{'='*20}")
        print(f" Running Example {ex_name} ")
        print(f"{'='*20}\n")
        EXAMPLES[ex_name]()

    plt.show()


if __name__ == "__main__":
    main()
