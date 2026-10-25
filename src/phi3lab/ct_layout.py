"""One saved geometry for HTML/SVG and editable TikZ, independent of identity."""
from collections import defaultdict
import hashlib,html,json,math
from .ct_identity import graph_id


def anchors(g):
    labels=[l for l,_ in g.external]
    if len(labels)==3:return dict(zip(labels,[(180,15),(20,260),(340,260)]))
    if len(labels)==2:return dict(zip(labels,[(15,145),(345,145)]))
    if len(labels)==1:return {labels[0]:(180,15)}
    return {l:(180+155*math.cos(-math.pi/2+2*math.pi*i/len(labels)),145+120*math.sin(-math.pi/2+2*math.pi*i/len(labels))) for i,l in enumerate(labels)}


def layout(g,mode='shared'):
    if g.n==1:return {0:(180,145)}
    if mode=='legacy-circle':
        return {i:(180+86*math.cos(-math.pi/2+2*math.pi*i/g.n),145+86*math.sin(-math.pi/2+2*math.pi*i/g.n)) for i in range(g.n)}
    if len(g.external)==3:fixed=dict(zip([l for l,_ in g.external],[(180,55),(65,220),(295,220)]))
    elif len(g.external)==2:fixed=dict(zip([l for l,_ in g.external],[(65,145),(295,145)]))
    else:fixed={l:(180+(x-180)*.65,145+(y-145)*.65) for l,(x,y) in anchors(g).items()}
    p={}
    for l,v in g.external:
        points=[fixed[k] for k,w in g.external if w==v]
        p[v]=tuple(sum(x[i] for x in points)/len(points) for i in (0,1))
    remaining=[v for v in range(g.n) if v not in p]
    for j,v in enumerate(remaining):
        if len(remaining)<=2 and len(g.external)==3:p[v]=(180+(j-(len(remaining)-1)/2)*90,140)
        else:p[v]=(180+65*math.cos(-math.pi/2+2*math.pi*j/max(1,len(remaining))),145+65*math.sin(-math.pi/2+2*math.pi*j/max(1,len(remaining))))
    for v in remaining:
        if g.vertex_kinds[v].startswith('C2') and len(g.external)==3:
            near=[b if a==v else a for a,b in g.edges if v in (a,b) and a!=b]
            if near:p[v]=tuple(sum(p[w][k] for w in near)/len(near) for k in (0,1))
    return p


def geometry(g,positions=None,external_positions=None,mode='shared',provenance=None):
    p=positions if positions is not None else layout(g,mode)
    if set(p)!=set(range(g.n)):raise ValueError('Layout requires every vertex exactly once')
    groups=defaultdict(list)
    for i,e in enumerate(g.edges):groups[e].append(i)
    curves={}
    for (a,b),ids in groups.items():
        x,y=p[a];xx,yy=p[b]
        for j,i in enumerate(ids):
            if a==b:
                radius=55+25*j
                points=[(x-4,y),(x-radius,y-radius),(x+radius,y-radius),(x+4,y)]
                label=(x,y-radius*.75)
            else:
                dx,dy=xx-x,yy-y;length=math.hypot(dx,dy) or 1;off=(j-(len(ids)-1)/2)*42
                # Do not draw a propagator through a third vertex after a
                # centroid contraction. That can visually erase an entire loop.
                if off==0 and length>0:
                    blocked=any(0.08<((vx-x)*dx+(vy-y)*dy)/(length*length)<0.92
                                and abs((vx-x)*dy-(vy-y)*dx)/length<17
                                for v,(vx,vy) in p.items() if v not in (a,b))
                    if blocked:off=70
                mx,my=(x+xx)/2-off*dy/length,(y+yy)/2+off*dx/length
                points=[(x,y),(x+2*(mx-x)/3,y+2*(my-y)/3),(xx+2*(mx-xx)/3,yy+2*(my-yy)/3),(xx,yy)]
                label=((x+xx)/2-off*dy/length/2+8*dy/length,(y+yy)/2+off*dx/length/2-8*dx/length)
            curves[str(i)]={'points':points,'label':label}
    rec={'schema':'phi3lab.shared-layout.v1','graph_id':graph_id(g),'mode':mode,'vertices':{str(v):p[v] for v in range(g.n)},
         'external':{str(l):xy for l,xy in (external_positions or anchors(g)).items()},'edges':curves,
         'provenance':provenance,'coordinates_are_display_only':True}
    rec['layout_id']=hashlib.sha256(json.dumps(rec,sort_keys=True).encode()).hexdigest()
    return rec


def svg(g,geo=None,marked=(),selected=None,bridges=()):
    geo=geo or geometry(g);p=geo['vertices']
    out=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 360 285" role="img" aria-label="φ³">',
         '<metadata>'+html.escape(json.dumps(geo))+'</metadata>',
         '<style>text{font:11px system-ui;fill:#163b45;text-anchor:middle}.edge{fill:none;stroke:#234b57;stroke-width:2.3}.marked{stroke:#d16d26;stroke-width:4}.bridge{stroke:#9e57a1;stroke-width:3}.selected{stroke:#c63232;stroke-width:3;stroke-dasharray:6 4}.el{font-size:10px;paint-order:stroke;stroke:white;stroke-width:4}</style>']
    for i,e in geo['edges'].items():
        ps=e['points'];x,y=e['label'];style='selected' if int(i)==selected else 'bridge' if int(i) in bridges else 'marked' if int(i) in marked else ''
        path='M'+','.join(map(str,ps[0]))+' C'+' '.join(','.join(map(str,xy)) for xy in ps[1:])
        out.append(f'<path class="edge {style}" data-eid="{i}" d="{path}"/><text class="el" x="{x}" y="{y}">e{i}</text>')
    for label,v in g.external:
        x,y=p[str(v)];xx,yy=geo['external'][str(label)]
        out.append(f'<path class="edge" d="M{x},{y} L{xx},{yy}"/><text x="{xx}" y="{yy-6}">p{label}</text>')
    for v in range(g.n):
        x,y=p[str(v)];kind=g.vertex_kinds[v]
        if kind!='phi3':out.append(f'<rect x="{x-16}" y="{y-11}" width="32" height="22" fill="#fff2cb" stroke="#926b1e"/><text x="{x}" y="{y+4}">{html.escape(kind)}</text>')
        else:out.append(f'<circle cx="{x}" cy="{y}" r="4" fill="#234b57"/><text x="{x+10}" y="{y+13}">v{v}</text>')
    return ''.join(out)+'</svg>'


def tikz(g,geo=None,marked=(),selected=None,bridges=()):
    geo=geo or geometry(g)
    coord=lambda xy:f'({xy[0]/55:.8f},{(285-xy[1])/55:.8f})'
    out=['% shared-layout-id: '+geo['layout_id'],r'\begin{tikzpicture}[scale=.65,every node/.style={font=\scriptsize}]']
    for v in range(g.n):
        ct=g.vertex_kinds[v]!='phi3';shape='rectangle,draw,fill=yellow!15' if ct else 'circle,fill,inner sep=1.5pt'
        out.append(r'\node['+shape+'] (v'+str(v)+') at '+coord(geo['vertices'][str(v)])+' {'+(g.vertex_kinds[v] if ct else '')+'};')
    for i,e in geo['edges'].items():
        color='red,dashed' if int(i)==selected else 'violet,thick' if int(i) in bridges else 'orange,thick' if int(i) in marked else 'black'
        ps=e['points'];out.append(r'\draw['+color+'] '+coord(ps[0])+' .. controls '+coord(ps[1])+' and '+coord(ps[2])+' .. '+coord(ps[3])+';')
        out.append(r'\node at '+coord(e['label'])+' {$e_{'+i+'}$};')
    for label,v in g.external:out.append(r'\draw (v'+str(v)+') -- '+coord(geo['external'][str(label)])+' node[above] {$p_{'+str(label)+'}$};')
    return '\n'.join(out+[r'\end{tikzpicture}'])


def export_family(rows,columns,out):
    """Export parent, marked subgraph, quotient and Henge from common geometry."""
    from pathlib import Path
    from .graphs import Graph
    root=Path(out);records={};names={c['id']:c.get('name',c['id']) for c in columns}
    (root/'figures').mkdir(exist_ok=True);(root/'latex').mkdir(exist_ok=True)
    def emit(name,g,p=None,marked=(),selected=None,bridges=(),provenance=None):
        geo=geometry(g,p,provenance=provenance);records[name]={'graph':g.to_dict(),'geometry':geo}
        (root/'figures'/f'{name}.svg').write_text(svg(g,geo,marked,selected,bridges),encoding='utf-8')
        (root/'latex'/f'{name}.tex').write_text(tikz(g,geo,marked,selected,bridges),encoding='utf-8')
    for row in rows:
        g=Graph.from_dict(row['graph']);name=row['old_id'];p=layout(g);emit(name,g,p)
        for h in row.get('henges',[]):
            suffix=f'{name}-henge-e{h["minimum_edge"]}';vm=h['original_vertex_to_skeleton']
            marked=[e for block in h['blocks'] if block['loops'] for e in block['edge_ids']]
            emit(suffix,g,p,marked,h['minimum_edge'],h['bridges'],{'parent':name,'edge_map':h['original_edge_map']})
            sg=Graph(len(h['blocks']),tuple(tuple(e) for e in h['skeleton_edges']),tuple((l,vm[v]) for l,v in g.external),tuple('H' if block['loops'] else 'phi3' for block in h['blocks']))
            pos={i:tuple(sum(p[v][k] for v in block['vertices'])/len(block['vertices']) for k in (0,1)) for i,block in enumerate(h['blocks'])}
            emit(f'{name}-skeleton-e{h["minimum_edge"]}',sg,pos,provenance={'parent':name,'vertices':vm,'edges':h['original_edge_map'],'placement':'centroid; retained external anchors'})
        for term in row.get('terms',[]):
            if term['kind']!='insertion' or 'host' not in term:continue
            label=names[term['column']];sub=next(s for s in row['recursive_ast']['proper_parts'] if s['edge_ids']==term['embedded_edges'])
            vertices=sub['vertices'];raw={0:tuple(sum(p[v][k] for v in vertices)/len(vertices) for k in (0,1))}
            keep=[v for v in range(g.n) if v not in vertices]
            for j,v in enumerate(keep,1):raw[j]=p[v]
            final=term['contraction_mapping']['vertices'];pos={final[v]:xy for v,xy in raw.items()}
            vm={v:final[0] for v in vertices};vm.update({v:final[j] for j,v in enumerate(keep,1)})
            provenance={'parent':name,'parent_to_host_vertices':vm,'contraction_mapping':term['contraction_mapping'],
                        'expanded_source_edges':term['embedded_edges'],'placement':'centroid; surviving vertices and external anchors unchanged'}
            emit(f'{name}-{label}-subgraph',g,p,term['embedded_edges'],provenance=provenance)
            emit(f'{name}-{label}-host',Graph.from_dict(term['host']),pos,provenance=provenance)
    (root/'layouts.json').write_text(json.dumps({'schema':'phi3lab.layout-bundle.v1','figures':records},indent=2),encoding='utf-8')
    return records
