"""Reproducible graph-to-executable-density artefacts through four loops."""
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from .graphs import generate as legacy_generate
from .higher_graphs import generate as higher_generate,symmetry
from .automatic_bph import compile_bph
from .provenance import source_fingerprint


def generated_graph(loops,legs,index=1):
    if type(loops)is not int or loops not in (1,2,3,4) or type(legs)is not int or legs not in (1,2,3):raise ValueError('Use1--4 loops and1--3 external legs')
    graphs=(legacy_generate if loops<=2 else higher_generate)(loops,legs)
    if type(index)is not int or not 1<=index<=len(graphs):raise ValueError('Graph index must lie between1 and'+str(len(graphs)))
    return graphs[index-1]


def expression_tex(compiled,ident='G'):
    """Executable coefficient plan written explicitly, including every edge row."""
    d=compiled.to_dict();lines=[r'\section*{'+ident+r': graph-derived recursive BPH density}',
        r'\[J_R=\int_{\alpha_e>0,\,\sum\alpha_e=1}\! d^{I-1}\alpha\;\mathcal K_G(\alpha).\]',
        r'\[\mathcal K_G=\sum_w(-1)^{|w|}\sum_{0\le k_\gamma\le\omega_\gamma}[t^k]U_w(t)^{-3}\Phi_E(H_w(t),M),\]',
        r'\[\Phi_3=-\log(1+H/M),\qquad\Phi_2=(M+H)\log(1+H/M)-H.\]',
        r'\[Q_w=\sum_e\alpha_e c_{e,w}c_{e,w}^{T},\quad U_w=\det Q_{\ell\ell},\quad H_w=\operatorname{Gram}(p):(Q_{pp}-Q_{p\ell}Q_{\ell\ell}^{-1}Q_{\ell p}).\]']
    if len(compiled.graph.external)==1:lines.append(r'\[J_R=0\quad\text{by the declared full one-point Taylor condition}.\]')
    for i,term in enumerate(d['terms']):
        lines.extend([r'\subsection*{Recursive word '+str(i)+r', sign '+str(term['sign'])+'}',
            r'Inner-first embedded parts: '+str(term['subgraphs_inner_first'])+r'. Jet orders: '+str(term['jet_orders'])+'.',r'\[\begin{aligned}'])
        for e,row in enumerate(term['scaled_routing']):
            entries=[]
            for v in row:
                c=v['coefficient']
                if not c:entries.append('0');continue
                monomial=''.join('t_{'+str(j+1)+'}'+('^{'+str(n)+'}' if n>1 else '') for j,n in enumerate(v['powers']) if n)
                entries.append(str(c)+monomial)
            lines.append('c_{'+str(e)+'}&=('+','.join(entries)+r')\\')
        lines.append(r'\end{aligned}\]')
    return '\n'.join(lines)+'\n'


def build_automatic(folder,language=None):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    (folder/'expressions').mkdir(exist_ok=True)
    fp=source_fingerprint();rows=[];groups=[];artifacts=[]
    from .automatic_bounds import primitive_uv_certificate
    for loops in range(1,5):
        for legs in (1,2,3):
            graphs=(legacy_generate if loops<=2 else higher_generate)(loops,legs)
            groups.append({'loops':loops,'legs':legs,'count':len(graphs)})
            for index,g in enumerate(graphs,1):
                ident=f'L{loops}E{legs}G{index}';compiled=compile_bph(g)
                program=compiled.to_dict();program.update(graph_id=ident,source_fingerprint=fp)
                certificate=primitive_uv_certificate(g) if legs==3 else {}
                if certificate.get('primitive'):program['primitive_uv_certificate']=certificate
                path=folder/'expressions'/(ident+'.json');path.write_text(json.dumps(program,ensure_ascii=False,indent=2))
                tex=path.with_suffix('.tex');tex.write_text(expression_tex(compiled,ident))
                alpha=[1/len(g.edges)]*len(g.edges)
                point=compiled.density(alpha)
                row={'id':ident,'loops':loops,'legs':legs,'index':index,'vertices':g.n,'edges':len(g.edges),
                     'symmetry_weight':'1' if symmetry(g)==1 else '1/'+str(symmetry(g)),
                     'proper_uv_parts':len(compiled.parts),'recursive_terms':len(compiled.terms),
                     'maximum_jet_dimension':max(map(lambda t:len(t['jet_orders']),compiled.terms),default=0),
                     'parameter_dimension':len(g.edges)-1,'equal_parameter_density':point,
                     'primitive':not bool(compiled.parts),'expression':'expressions/'+path.name,
                     'tex':'expressions/'+tex.name,'expression_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                     'graph':g.to_dict()}
                if certificate.get('primitive'):row['rigorous_primitive_bound_coefficient']=certificate.get('hepp_bound')
                rows.append(row);artifacts.extend((path,tex))
    data={'schema_version':1,'generated_utc':datetime.now(timezone.utc).isoformat(),'source_fingerprint':fp,
          'source_unchanged':fp==source_fingerprint(),'method':'Recursive BPH graph compiler; no topology-specific integral dispatch',
          'groups':groups,'total_graphs':len(rows),'compiled_graphs':len(rows),'graphs':rows,
          'all_ok':len(rows)==698 and fp==source_fingerprint(),
          'scope':'All1--4-loop fixed-external-label connected amputated trivalent 1PI graphs with E=1,2,3 in massive Euclidean6D phi3; all expressions compiled and equal-parameter points evaluated. Integral execution coverage is separately recorded.'}
    path=folder/'catalog.json';path.write_text(json.dumps(data,ensure_ascii=False,indent=2));artifacts.append(path)
    locale_path=Path(__file__).resolve().parents[1]/'release-language.json'
    if language is None and locale_path.exists():language=json.loads(locale_path.read_text())['language']
    page=folder/'index.html'
    if language in ('cn','en'):
        from .localized_reporting import write_localized_html,localize_expression_tree
        localize_expression_tree(folder/'expressions',language)
        publication_path=folder.parent/'release-summary.json'
        publication=json.loads(publication_path.read_text()) if publication_path.exists() else {}
        write_localized_html(data,page,language=language,publication=publication)
    else:
        from .automatic_reporting import write_automatic_html
        write_automatic_html(data,page)
    artifacts.append(page)
    manifest={'source_fingerprint':fp,'files':{str(p.relative_to(folder)):hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts}}
    (folder/'manifest.json').write_text(json.dumps(manifest,indent=2))
    return data


def evaluate_automatic(loops=2,legs=3,index=1,*,graph=None,invariants=None,m2=1,integrate=False,power=14,replicates=8,beta=.35,sampler='auto'):
    g=graph if graph is not None else generated_graph(loops,legs,index)
    compiled=compile_bph(g);ident='input' if graph is not None else f'L{loops}E{legs}G{index}'
    result={'graph_id':ident,'source_fingerprint':source_fingerprint(),'expression':compiled.to_dict(),
            'equal_parameter_density':compiled.density([1/len(g.edges)]*len(g.edges),invariants,m2),
            'invariants':invariants,'m2':m2}
    if integrate:
        from .automatic_integrate import integrate_compiled
        result['integral']=integrate_compiled(compiled,invariants,m2,power=power,replicates=replicates,beta=beta,sampler=sampler)
    return result
