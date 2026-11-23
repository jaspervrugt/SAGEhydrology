"""Build isolated CaRT/CaRT+ continuous solver libraries and provenance."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess


def build():
    root = Path(__file__).resolve().parent
    compiler = shutil.which('g++') or r'C:\Strawberry\c\bin\g++.exe'
    outputs = {}
    for folder, source, name, definitions in (
        ('cart', 'cart_ode.cpp', 'cart_numerics_v1.dll', []),
        ('CaRTPlus', 'cartplus_continuous.cpp', 'cartplus_numerics_v1.dll', []),
        ('CaRTPlus', 'cartplus_continuous.cpp', 'cartplus_numerics_v1_35.dll',
         ['-DCARTPLUS_CALIBRATED_ONLY']),
    ):
        target = root / folder / name
        subprocess.run([compiler, '-std=c++17', '-O3', '-shared', '-static',
                        '-DCART_FAST_SCALAR', *definitions,
                        str(root / folder / source), '-o', str(target)], check=True)
        outputs[str(target.relative_to(root))] = hashlib.sha256(target.read_bytes()).hexdigest()
    sources = [root/'cart/cart_ode.cpp', root/'CaRTPlus/cartplus_continuous.cpp',
               *sorted((root/'CaRTPlus').glob('*.hpp'))]
    manifest = dict(revision='cart-numerics-v1', libraries=outputs,
                    sources={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sources})
    (root/'cart_numerics_build.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(outputs, indent=2))


if __name__ == '__main__':
    build()
