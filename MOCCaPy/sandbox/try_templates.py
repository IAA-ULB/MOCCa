"""
Experiment with string.Template
This might be useful for in Hephaestos and its offspring.

+ templates are basically strings with variables that can be replaced.
+ everything is Python
+ no template files that can be mistaken for Python files

- no syntax highlighting when editing the templates (since they are just Python strings)
- Code itself becomes harder to read, boundaries between template and code is less clear

The lack of syntax highlighting is problematic in view of the size of the templates.
"""
from string import Template
template = Template('''
"""
    Relativity. Created by `try_templates.py`
"""

lines = [
    "There was a young lady named $name",
    "Whose ${property} was far faster than light",
    "She set out one day",
    "In a relative way"
    "And returned on the previous night.",
]
print('\\n'.join(lines))
''')

if __name__ == "__main__":
    with open("relativity.py", "w") as f:
        print(template.substitute(name='Bright', property='speed'), file=f)

    import relativity # prints the limerick, just to make sure everything works as intended.
