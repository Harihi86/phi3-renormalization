"""Explicit, optional pySecDec bridge. The standard-library core stays independent.

The exports are BARE integrals in d=D-2*eps, not integrated forest sums.
All kinematic invariants s_i are Euclidean p_i**2; pySecDec uses Minkowski
scalar products, hence the negative squares in the replacement rules.
"""
import json
from pathlib import Path


def pysecdec_spec(g, graph_id, dimension=6):
    if dimension not in (4, 6):
        raise ValueError("Use dimension 4 or 6")
    labels = [label for label, _ in g.external]
    if not 1 <= len(labels) <= 3:
        raise ValueError("The bridge supports one to three external legs")
    momenta = [f'p{x}' for x in labels]
    rules = [(momenta[-1], '-'+ '-'.join(momenta[:-1]) if len(labels)>1 else '0')]
    parameters = ['m2']
    if len(labels) >= 2:
        parameters.append('s1')
        rules.append((momenta[0]+'*'+momenta[0], '-s1'))
    if len(labels) == 3:
        parameters.extend(['s2', 's3'])
        rules.extend([(momenta[1]+'*'+momenta[1], '-s2'),
                      (momenta[0]+'*'+momenta[1], '(s1+s2-s3)/2')])
    rules.append(('m**2', 'm2'))
    name = 'phi3_'+''.join(c.lower() if c.isalnum() else '_' for c in graph_id)
    return {
        'graph_id': graph_id, 'package_name': name,
        'loop_integral': {
            'internal_lines': [['m', list(e)] for e in g.edges],
            'external_lines': [[f'p{label}', v] for label, v in g.external],
            'replacement_rules': [list(r) for r in rules],
            'Feynman_parameters': [f'x{i}' for i in range(len(g.edges))],
            'dimensionality': f'{dimension}-2*eps',
        },
        'real_parameters': parameters,
        'example_point': [1]*len(parameters),
        'requested_orders': [0],
        'additional_prefactor': (-1)**len(g.edges),
        'normalization': {
            'target': 'J_E = integral prod_l[d^d k_l / pi^(d/2)] prod_e[1/(q_e^2+m2)]',
            'mass': 'm2 is the positive mass SQUARED; every internal line has mass sqrt(m2).',
            'kinematics': 's_i=p_i,E^2; p_last=-sum(other external momenta). For E=2 only s1 is independent.',
            'minkowski_conversion': 'p_i,M dot p_j,M = -p_i,E dot p_j,E; additional_prefactor=(-1)^I cancels the pySecDec unit-power overall sign.',
            'excluded': ['couplings', 'diagram symmetry weight 1/S', '(4*pi)^(-L*d/2)', 'renormalization-scale factors', 'BPH forest subtraction', 'MS or MSbar finite factors'],
        },
        'status': 'Runnable input export; generation, compilation and integration are separate operations. No integral value is asserted by export.',
    }


_GENERATE = '''"""Inspect by default; --generate writes C++ sector-decomposition code."""
import argparse,json,os
from pathlib import Path
import pySecDec as psd

root=Path(__file__).resolve().parent
config=json.loads((root/'integral.json').read_text())
p=argparse.ArgumentParser(description=config['status'])
p.add_argument('--generate',action='store_true')
args=p.parse_args()
li=psd.LoopIntegralFromGraph(**config['loop_integral'])
print('pySecDec',psd.__version__)
print('U =',li.U)
print('F =',li.F)
print('Euclidean prefactor =',config['additional_prefactor']*li.Gamma_factor)
print('This is a BARE dimensionally regulated integral, without 1/S or couplings.')
if args.generate:
    build=root/'build';build.mkdir(exist_ok=True);os.chdir(build)
    psd.loop_package(name=config['package_name'],loop_integral=li,
        requested_orders=config['requested_orders'],
        real_parameters=config['real_parameters'],
        additional_prefactor=config['additional_prefactor'],
        contour_deformation=False)
'''

_INTEGRATE = '''"""Run only after generate.py --generate and compiling the emitted package."""
import json
from pathlib import Path
from pySecDec.integral_interface import IntegralLibrary

root=Path(__file__).resolve().parent
c=json.loads((root/'integral.json').read_text())
library=root/'build'/c['package_name']/(c['package_name']+'_pylink.so')
if not library.exists():
    raise SystemExit('Compile the generated library first; see ../README.txt')
li=IntegralLibrary(str(library));li.use_Qmc()
result=li(real_parameters=c['example_point'],epsrel=1e-5,epsabs=1e-8)
print('parameters:',dict(zip(c['real_parameters'],c['example_point'])))
print('without prefactor / prefactor / with prefactor:',*result,sep='\\n')
print('Bare J_E Laurent series: not a BPH-renormalized amplitude.')
'''


def write_pysecdec_export(data, folder):
    from .graphs import Graph
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for record in data['graphs']:
        target = folder/record['id']
        target.mkdir(exist_ok=True)
        spec = pysecdec_spec(Graph.from_dict(record), record['id'], data['dimension'])
        spec['diagram_weight_excluded'] = record['weight']
        for name, content in [('integral.json', json.dumps(spec, ensure_ascii=False, indent=2)),
                              ('generate.py', _GENERATE), ('integrate.py', _INTEGRATE)]:
            path = target/name
            path.write_text(content, encoding='utf-8')
            written.append(path)
    readme = folder/'README.txt'
    readme.write_text('''Phi3Lab → pySecDec 1.6.6\n\n每张图的 integral.json 都是明确约定的裸积分输入；核心无需安装 pySecDec。\n1. 在独立环境安装官方 pySecDec 1.6.6。Linux 有官方预编译包；Apple Silicon 可能需要源码构建及其依赖。\n2. 进入某张图目录，运行 python generate.py，先核查 U、F、Gamma 因子。\n3. 运行 python generate.py --generate；需要完整 pySecDec 及配套工具。\n4. 在 build/phi3_<图编号小写> 中运行 GNU make -j2。官方工具链要求 C++17。\n5. 回到图目录运行 python integrate.py，得到含误差的裸积分 Laurent 系列。\n\n所有 s_i 为 Euclidean 外腿平方；默认 m2=s_i=1 是合法 Euclidean 点。\n目标测度 d^d k/pi^(d/2)，d=D−2ε。additional_prefactor=(-1)^I 已补偿 pySecDec 的 Minkowski 整体符号。\n图权重、耦合、(4π)^(-Ld/2)、μ 因子、MSbar 因子均未乘入；本输出没有执行 BPH 减除。\n不要将裸 ε 极点与有质量零动量减除后的有限结果直接比较。\n自动互校工具只需上游 Python 源码与 Python 依赖，核查 U/F、维数和 Gamma 因子；这不等于 C++ 生成、编译或积分已完成。\n官方文档：https://secdec.readthedocs.io/en/stable/gettingstarted.html\n源码许可证 GPL-3.0；本项目未打包该源码。\n''', encoding='utf-8')
    written.append(readme)
    return written
