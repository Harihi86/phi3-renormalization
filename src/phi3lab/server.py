"""Loopback-only standard-library workbench. No network service is required."""
import json,threading,webbrowser
from functools import partial
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs,urlsplit
from .workbench import evaluate_overlap
from .two_loop import evaluate_two_loop
from .two_loop_registry import two_loop_catalog


class WorkbenchHandler(SimpleHTTPRequestHandler):
    # Limit simultaneous high-order cubature requests from the browser.
    calculation_slots=threading.BoundedSemaphore(2)

    def _json(self,status,value):
        body=json.dumps(value,ensure_ascii=False,allow_nan=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type','application/json; charset=utf-8')
        self.send_header('Content-Length',str(len(body)))
        self.send_header('Cache-Control','no-store')
        self.end_headers();self.wfile.write(body)

    def do_GET(self):
        parsed=urlsplit(self.path)
        if parsed.path=='/api/health':
            self._json(200,{'service':'Phi3Lab BPH workbench','ok':True});return
        if parsed.path not in ('/api/overlap','/api/two-loop','/api/automatic'):
            if parsed.path.startswith('/api/'):
                self._json(404,{'error':'未知计算入口'});return
            return super().do_GET()
        if not self.calculation_slots.acquire(blocking=False):
            self._json(503,{'error':'已有两个计算正在进行，请稍后重试'});return
        try:
            parameters=parse_qs(parsed.query,keep_blank_values=True)
            if any(len(v)!=1 for v in parameters.values()):raise ValueError('每个参数只能填写一次')
            allowed=({'s','m2','order'} if parsed.path=='/api/overlap' else
                     {'loops','legs','index','s1','s2','s3','m2','integrate','power','replicates','beta','sampler'} if parsed.path=='/api/automatic' else
                     {'graph','s1','s2','s3','m2','order'})
            unknown=set(parameters)-allowed
            if unknown:raise ValueError('未知参数：'+', '.join(sorted(unknown)))
            m2=float(parameters.get('m2',['1'])[0])
            if parsed.path=='/api/automatic':
                from .automatic_catalog import evaluate_automatic
                loops=int(parameters.get('loops',['2'])[0]);legs=int(parameters.get('legs',['3'])[0]);index=int(parameters.get('index',['1'])[0])
                n=0 if legs==1 else 1 if legs==2 else 3
                if any('s'+str(j) in parameters for j in range(n+1,4)):raise ValueError('此图不需要这些额外的不变量')
                invariants=[float(parameters.get('s'+str(j),['1'])[0]) for j in range(1,n+1)]
                integrate=parameters.get('integrate',['0'])[0]
                if integrate not in ('0','1'):raise ValueError('integrate须为0或1')
                power=int(parameters.get('power',['12'])[0]);replicates=int(parameters.get('replicates',['4'])[0])
                if not 6<=power<=18 or not 2<=replicates<=16:raise ValueError('页面求积范围为2^6--2^18点、2--16组；更大运行请使用命令行')
                beta=float(parameters.get('beta',['.35'])[0]);sampler=parameters.get('sampler',['auto'])[0]
                result=evaluate_automatic(loops,legs,index,invariants=invariants,m2=m2,integrate=integrate=='1',power=power,replicates=replicates,beta=beta,sampler=sampler)
            elif parsed.path=='/api/overlap':
                s=float(parameters.get('s',['1'])[0])
                order=int(parameters.get('order',['16'])[0])
                result=evaluate_overlap(s,m2,order)
            else:
                graph=parameters.get('graph',['L2E2G1'])[0]
                catalog={r['id']:r for r in two_loop_catalog()}
                if graph not in catalog:raise ValueError('请选择有效的两圈图编号')
                legs=len(catalog[graph]['external_labels'])
                n=0 if legs==1 else 1 if legs==2 else 3
                if any('s'+str(j) in parameters for j in range(n+1,4)):
                    raise ValueError('此图不需要这些额外的动量不变量')
                invariants=[float(parameters.get('s'+str(j),['1'])[0]) for j in range(1,n+1)]
                order=int(parameters['order'][0]) if 'order' in parameters else None
                result=evaluate_two_loop(graph,invariants,m2,order)
            self._json(200,result)
        except (ValueError,OverflowError,ArithmeticError,RuntimeError) as exc:
            self._json(400,{'error':str(exc)})
        finally:
            self.calculation_slots.release()

    def log_message(self,format,*args):
        # A local scientific workbench needs no persistent access log.
        pass


def make_server(directory,port=8767):
    if type(port) is not int or not 0 <= port <= 65535:raise ValueError('端口须为 0 到 65535 的整数')
    root=Path(directory).resolve()
    if not (root/'results/index.html').is_file():raise ValueError('未找到结果页，请先运行 build --wick')
    return ThreadingHTTPServer(('127.0.0.1',port),partial(WorkbenchHandler,directory=str(root)))


def serve(directory,port=8767,open_browser=True):
    server=make_server(directory,port)
    url=f'http://127.0.0.1:{server.server_port}/results/index.html#two-loop-workbench'
    print('本地 BPH 工作台：'+url,flush=True)
    print('保持此窗口开启；按 Ctrl+C 结束。计算在本机进行。',flush=True)
    if open_browser:webbrowser.open(url)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
