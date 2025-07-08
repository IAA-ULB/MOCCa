#!/usr/bin/env python
# a simple python script to look for link libs.

import os
import click


@click.command()
@click.option('-f', '--find', default="", help='list all entries containing name.')
@click.option('--so', is_flag=True, help='List only files ending on ".so"')
def scan(find, so):
    # create a list of directories in LD_LIBRARY_PATH, where the linker will
    # look for needed libraries.
    libs = os.environ["LD_LIBRARY_PATH"].split(':')
    print(f"[[{find}]]")
    found = False
    for d in libs:
        entries = os.listdir(d)
        for e in entries:
            # print(f"{d}: {e}")
            if find in e:
                if not so or e.endswith('so'):
                    print(f"{d}: {e}")
                    found = True

    if not found:
        print("Nothing found.")


if __name__ == "__main__":
    scan()