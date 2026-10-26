'use strict';
(() => {
const $=id=>document.getElementById(id), W=900,H=590, TMAX=300;
const config={pore:300,curvature:-0.07,stiffness:1,adhesion:.46,motility:.68,division:.56,geometry:'lattice',lineage:'osteogenic',chem:false,count:65,seed:42};
const sliders=[['pore','Pore diameter',100,600,10,'µm',0],['curvature','Surface curvature',-.2,.2,.01,'µm⁻¹',2],['stiffness','Stiffness',.5,2,.1,'GPa',1],['adhesion','Cell adhesion',0,1,.02,'relative',2],['motility','Cell motility',0,1,.02,'relative',2],['division','Proliferation propensity',0,1,.02,'relative',2]];
let thresholds={coverage:30,commitment:20,infiltration:38};
let exp=null, history=[],running=false,lastNow=0,acc=0;
const sctx=$('simCanvas').getContext('2d'), cctx=$('chartCanvas').getContext('2d'), mctx=$('morphCanvas').getContext('2d');
const presets={colonization:{desc:'How does pore architecture influence migration and infiltration?',lineage:'vascular',pore:390,curvature:0.02},fate:{desc:'How do mechanical and geometric cues bias cell differentiation?',lineage:'osteogenic',pore:280,curvature:-.08},kernel:{desc:'Run a parameter sweep and classify biologically acceptable configurations.',lineage:'osteogenic',pore:280,curvature:-.07}};
function rngFromSeed(seed){let s=(seed|0)||1;return()=>{s^=s<<13;s^=s>>>17;s^=s<<5;return (s>>>0)/4294967296;}}
function gauss(rng){return (rng()+rng()+rng()+rng()-2)*1.732;}
function clamp(v,a,b){return Math.max(a,Math.min(b,v));}
function round(v,n=1){return +v.toFixed(n);}
function setupControls(){
 const root=$('sliders');root.innerHTML=sliders.map(([key,name,min,max,step,unit,n])=>`<div class="slider-row"><div class="label-line"><label for="s_${key}">${name}</label><output id="o_${key}">${(+config[key]).toFixed(n)} ${unit}</output></div><input type="range" id="s_${key}" min="${min}" max="${max}" step="${step}" value="${config[key]}"></div>`).join('');
 sliders.forEach(([key,,_,__,___,unit,n])=>$('s_'+key).addEventListener('input',()=>{config[key]=+$('s_'+key).value;$('o_'+key).textContent=config[key].toFixed(n)+' '+unit;drawMorph();}));
 $('thresholds').innerHTML=[['coverage','Spatial coverage',10,80,'%'],['commitment','Committed fraction',0,85,'%'],['infiltration','Infiltration',10,90,'%']].map(([key,name,min,max,u])=>`<div class="slider-row"><div class="label-line"><label>${name}</label><output id="t_${key}">${thresholds[key]}${u}</output></div><input type="range" min="${min}" max="${max}" value="${thresholds[key]}" id="r_${key}"></div>`).join('');
 Object.keys(thresholds).forEach(k=>$('r_'+k).addEventListener('input',()=>{thresholds[k]=+$('r_'+k).value;$('t_'+k).textContent=thresholds[k]+'%';renderHistory();drawMorph();}));
 }
function readControls(){config.geometry=$('architecture').value;config.lineage=$('lineage').value;config.chem=$('chemEnabled').checked;config.count=+$('cellCount').value;config.seed=+$('seed').value||1;}
function writeControls(){sliders.forEach(([key,,a,b,c,unit,n])=>{$('s_'+key).value=config[key];$('o_'+key).textContent=(+config[key]).toFixed(n)+' '+unit;});$('architecture').value=config.geometry;$('lineage').value=config.lineage;$('chemEnabled').checked=config.chem;$('cellCount').value=String(config.count);$('seed').value=config.seed;}
function scaffoldSolid(x,y,c){
 const pitch=108;const narrow=clamp((c.pore-100)/500,0,1);const width=clamp(23-15*narrow,6,26);
 if(c.geometry==='channels')return Math.abs(((y+pitch/2)%pitch)-pitch/2)<width*1.4 && Math.sin(x/60)<.6;
 let xx=x,yy=y;if(c.geometry==='hex')xx+=Math.floor(y/(pitch*.9))%2*pitch*.5;
 if(c.geometry==='irregular'){xx+=15*Math.sin(y*.026+c.seed*.1);yy+=14*Math.sin(x*.028+c.seed*.04);}
 const dx=Math.abs(((xx+pitch/2)%pitch)-pitch/2),dy=Math.abs(((yy+pitch/2)%pitch)-pitch/2);
 // Scaffold struts at the vertices leave connected free paths for cells to navigate.
 return dx<width && dy<width;
}
function localGeo(x,y,c){const poreScore=clamp((c.pore-100)/500,0,1);let curv=1-Math.abs(c.curvature+0.07)/0.28;return clamp(.33+.38*poreScore+.18*curv+.13*Math.sin(x*.019+y*.023),.07,1.0);}
function makeExperiment(c){const rr=rngFromSeed(c.seed);const state={c:{...c},rr,t:0,cells:[],matrix:[],series:[],last:null};for(let i=0;i<c.count;i++){
  let x=28+rr()*110,y=25+rr()*(H-50),trial=0;while(scaffoldSolid(x,y,c)&&trial++<60){x=28+rr()*120;y=25+rr()*(H-50);}
  state.cells.push({x,y,age:rr()*80,phen:0,stress:0,path:[],phase:rr()*6.28});}
 snapshot(state);return state;}
function step(e){if(e.t>=TMAX)return;let c=e.c,rr=e.rr,old=e.cells.slice(),grid=new Map(),size=32;
 for(const cell of old){const gx=Math.floor(cell.x/size),gy=Math.floor(cell.y/size),k=gx+','+gy;if(!grid.has(k))grid.set(k,[]);grid.get(k).push(cell);}
 const babies=[];
 for(const cell of old){const gx=Math.floor(cell.x/size),gy=Math.floor(cell.y/size);let neighbors=[],fx=0,fy=0;
 for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++){for(const other of grid.get((gx+dx)+','+(gy+dy))||[]){if(other===cell)continue;const vx=cell.x-other.x,vy=cell.y-other.y,d=Math.hypot(vx,vy);if(d<28){neighbors.push(other);if(d<14&&d>.01){fx+=vx/d*(14-d)/14*.65;fy+=vy/d*(14-d)/14*.65;}if(d>13&&d<26){fx-=vx/d*c.adhesion*.04;fy-=vy/d*c.adhesion*.04;}}}}
 const crowd=clamp(neighbors.length/12,0,1);
 const g=localGeo(cell.x,cell.y,c),chem=c.chem?.4:0;
 const direction=rr()*Math.PI*2;
 let vx=(Math.cos(direction)*.62+fx+(.22+chem)*g)*(.4+c.motility*1.9);
 let vy=(Math.sin(direction)*.77+fy+Math.sin(cell.y*.016+cell.x*.009)*.12)*(.4+c.motility*1.7);
 if(c.geometry==='channels')vy*=.66;
 let nx=clamp(cell.x+vx,9,W-9),ny=clamp(cell.y+vy,9,H-9);
 if(scaffoldSolid(nx,ny,c)){if(!scaffoldSolid(nx,cell.y,c)){ny=cell.y;}else if(!scaffoldSolid(cell.x,ny,c)){nx=cell.x;}else{nx=cell.x;ny=cell.y;}}
 if(e.t%4===0){cell.path.push([cell.x,cell.y]);if(cell.path.length>11)cell.path.shift();}
 cell.x=nx;cell.y=ny;cell.age++;
 const strain=Math.abs(c.stiffness-({osteogenic:1.35,proliferative:.9,vascular:.65}[c.lineage]||1))/1.5;
 cell.stress=clamp(crowd*.62+strain*.28+(1-g)*.23,0,1);
 const targetPore={osteogenic:290,proliferative:200,vascular:470}[c.lineage];
 const geomFit=Math.exp(-Math.pow((c.pore-targetPore)/235,2)-Math.pow((c.curvature-({osteogenic:-.09,proliferative:0,vascular:.07}[c.lineage]))/.19,2));
 const stiffFit=Math.exp(-Math.pow((c.stiffness-({osteogenic:1.35,proliferative:.9,vascular:.7}[c.lineage]))/.9,2));
 const fateDrive=.20+geomFit*.48+stiffFit*.28;
 if(!cell.phen&&cell.age>26&&rr()<.0035*fateDrive*(1-cell.stress*.6))cell.phen=1;
 if(e.t%2===0 && cell.phen && rr()<.16){e.matrix.push({x:cell.x+(rr()-.5)*12,y:cell.y+(rr()-.5)*12});if(e.matrix.length>1600)e.matrix.shift();}
 const probability=.0042*c.division*(1-crowd*.9)*(1-cell.stress*.75)*(cell.phen?.36:1)*(1-.35*Math.abs(c.stiffness-1));
 if(cell.age>18&&old.length+babies.length<450&&rr()<probability){let bx=clamp(cell.x+gauss(rr)*14,8,W-8),by=clamp(cell.y+gauss(rr)*14,8,H-8);if(!scaffoldSolid(bx,by,c)){babies.push({x:bx,y:by,age:0,phen:0,stress:0,path:[],phase:rr()*6.28});cell.age=0;}}
 }
 e.cells.push(...babies);e.t++;if(e.t%4===0||e.t===TMAX)snapshot(e);
}
function outcomes(e){const living=e.cells.length,committed=e.cells.filter(c=>c.phen).length;const bins=new Set();let xmax=0,stress=0;
 for(const c of e.cells){bins.add(Math.floor(c.x/45)+','+Math.floor(c.y/45));xmax+=c.x;stress+=c.stress;}
 return {population:living,coverage:clamp(100*bins.size/260,0,100),commitment:living?100*committed/living:0,infiltration:living?100*(xmax/living)/W:0,matrix:Math.min(100,100*e.matrix.length/500),stress:living?100*stress/living:0};}
function snapshot(e){e.last=outcomes(e);e.series.push({t:e.t,...e.last});}
function selectedStatus(o){return o.coverage>=thresholds.coverage&&o.commitment>=thresholds.commitment&&o.infiltration>=thresholds.infiltration;}
function renderSimulation(){const ctx=sctx,c=exp.c;ctx.clearRect(0,0,W,H);ctx.fillStyle='#f5f9fa';ctx.fillRect(0,0,W,H);
 if($('showField').checked){for(let y=0;y<H;y+=24)for(let x=0;x<W;x+=24){const z=localGeo(x,y,c);ctx.fillStyle=`rgba(27,143,153,${z*.13})`;ctx.fillRect(x,y,24,24);}}
 ctx.fillStyle='#a6b5bb';for(let y=0;y<H;y+=3)for(let x=0;x<W;x+=3)if(scaffoldSolid(x,y,c))ctx.fillRect(x,y,3,3);
 ctx.fillStyle='rgba(144,105,187,.19)';for(const z of exp.matrix){ctx.beginPath();ctx.arc(z.x,z.y,3.2,0,Math.PI*2);ctx.fill();}
 for(const cell of exp.cells){if($('showPaths').checked&&cell.path.length){ctx.beginPath();ctx.moveTo(cell.path[0][0],cell.path[0][1]);for(const pt of cell.path)ctx.lineTo(...pt);ctx.lineTo(cell.x,cell.y);ctx.strokeStyle='rgba(36,107,151,.24)';ctx.lineWidth=1;ctx.stroke();}
 ctx.beginPath();ctx.arc(cell.x,cell.y,5.1,0,Math.PI*2);ctx.fillStyle=cell.stress>.78?'#e48e48':cell.phen?'#0b9c78':'#2f74aa';ctx.fill();ctx.lineWidth=1;ctx.strokeStyle='#ffffffb0';ctx.stroke();}
 ctx.fillStyle='#223e50';ctx.font='13px sans-serif';ctx.fillText('Cell seeding zone',17,22);ctx.strokeStyle='#627d8c77';ctx.setLineDash([5,5]);ctx.beginPath();ctx.moveTo(150,0);ctx.lineTo(150,H);ctx.stroke();ctx.setLineDash([]);
 $('timeLabel').textContent=`t = ${exp.t} / ${TMAX}`;
}
function renderMetrics(){const o=exp.last;const entries=[['Cells',o.population,''],['Coverage',round(o.coverage), '%'],['Commitment',round(o.commitment),'%'],['Infiltration',round(o.infiltration),'%'],['Matrix index',round(o.matrix),'%']];$('metrics').innerHTML=entries.map(([name,val,u])=>`<div class="metric"><span class="value">${val}${u}</span><span class="name">${name}</span></div>`).join('');$('statusPill').textContent=running?'Running':exp.t>=TMAX?'Complete':'Paused';$('viabilityText').textContent=`Current run: ${selectedStatus(o)?'meets':'does not meet'} the selected teaching thresholds. Coverage ${round(o.coverage)}%, commitment ${round(o.commitment)}%, infiltration ${round(o.infiltration)}%.`;
}
function graphBase(ctx,width,height,title,xlabel,ylabel){ctx.clearRect(0,0,width,height);ctx.fillStyle='#fff';ctx.fillRect(0,0,width,height);ctx.strokeStyle='#e5edf0';ctx.lineWidth=1;ctx.font='12px sans-serif';ctx.fillStyle='#627987';for(let i=0;i<=4;i++){const y=height-32-i*(height-54)/4;ctx.beginPath();ctx.moveTo(46,y);ctx.lineTo(width-15,y);ctx.stroke();ctx.fillText(String(i*25),10,y+4);}ctx.fillText(xlabel,width/2-25,height-7);if(title)ctx.fillText(title,55,15);}
function renderChart(){const ctx=cctx,w=900,h=235;graphBase(ctx,w,h,'Normalized outcomes over time','Time','');const series=exp.series;const entries=[['coverage','#248a96'],['commitment','#149e76'],['infiltration','#c87942']];for(const [key,color] of entries){ctx.beginPath();series.forEach((d,i)=>{const x=46+d.t/TMAX*(w-65),y=h-32-clamp(d[key],0,100)/100*(h-54);i?ctx.lineTo(x,y):ctx.moveTo(x,y);});ctx.strokeStyle=color;ctx.lineWidth=2.5;ctx.stroke();}ctx.font='12px sans-serif';let x=610;for(const [key,color] of entries){ctx.fillStyle=color;ctx.fillRect(x,11,12,3);ctx.fillText(key,x+16,16);x+=key==='commitment'?112:93;}}
function drawMorph(){const ctx=mctx,w=600,h=390;ctx.clearRect(0,0,w,h);ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);const x0=57,x1=580,y0=25,y1=340;
 const xOf=p=>x0+(p-100)/500*(x1-x0);const yOf=c=>y1-(c+.2)/.4*(y1-y0);
 // Clearly labelled reference region, not fitted to data
 ctx.fillStyle='rgba(54,158,134,.12)';ctx.strokeStyle='rgba(34,135,118,.38)';ctx.lineWidth=2;ctx.beginPath();ctx.ellipse(xOf(290),yOf(-.07),155,102,-.1,0,Math.PI*2);ctx.fill();ctx.stroke();
 ctx.strokeStyle='#e5edf0';ctx.lineWidth=1;ctx.font='12px sans-serif';ctx.fillStyle='#526b79';for(let p=100;p<=600;p+=100){let x=xOf(p);ctx.beginPath();ctx.moveTo(x,y0);ctx.lineTo(x,y1);ctx.stroke();ctx.fillText(String(p),x-11,y1+18);}for(let z=-.2;z<=.201;z+=.1){let y=yOf(z);ctx.beginPath();ctx.moveTo(x0,y);ctx.lineTo(x1,y);ctx.stroke();ctx.fillText(z.toFixed(1),15,y+4);}
 for(const item of history){ctx.beginPath();ctx.arc(xOf(item.c.pore),yOf(item.c.curvature),5.5,0,Math.PI*2);ctx.fillStyle=selectedStatus(item.o)?'#149978':'#c87942';ctx.fill();ctx.strokeStyle='#fff';ctx.lineWidth=1;ctx.stroke();}
 ctx.beginPath();ctx.arc(xOf(config.pore),yOf(config.curvature),9,0,Math.PI*2);ctx.fillStyle='#1c4777';ctx.fill();ctx.lineWidth=2;ctx.strokeStyle='white';ctx.stroke();ctx.fillStyle='#526b79';ctx.font='13px sans-serif';ctx.fillText('Pore diameter (µm)',240,382);ctx.save();ctx.translate(13,244);ctx.rotate(-Math.PI/2);ctx.fillText('Surface curvature (µm⁻¹)',0,0);ctx.restore();ctx.font='11px sans-serif';ctx.fillStyle='#137e7b';ctx.fillText('Hypothetical reference region',175,30);
}
function runToEnd(e,limit=TMAX){while(e.t<limit)step(e);return e.last;}
function makeRecord(e){return {id:history.length+1,date:new Date().toISOString(),c:{...e.c},t:e.t,o:{...e.last},viable:selectedStatus(e.last)};}
function renderHistory(){$('runCount').textContent=history.length+' saved';const root=$('runTable');if(!history.length){root.innerHTML='<p class="hint">No experiments saved yet.</p>';return;}
 root.innerHTML='<table><thead><tr><th>#</th><th>Architecture</th><th>Pore</th><th>Curvature</th><th>Lineage</th><th>Cells</th><th>Coverage</th><th>Committed</th><th>Infiltration</th><th>Viable?</th><th></th></tr></thead><tbody>'+history.map((r,i)=>`<tr><td>${i+1}</td><td>${r.c.geometry}</td><td>${r.c.pore}</td><td>${r.c.curvature.toFixed(2)}</td><td>${r.c.lineage}</td><td>${r.o.population}</td><td>${round(r.o.coverage)}%</td><td>${round(r.o.commitment)}%</td><td>${round(r.o.infiltration)}%</td><td>${selectedStatus(r.o)?'Yes':'No'}</td><td><button data-load="${i}">Load</button></td></tr>`).join('')+'</tbody></table>';
 root.querySelectorAll('[data-load]').forEach(b=>b.onclick=()=>{running=false;Object.assign(config,history[+b.dataset.load].c);writeControls();restart();});
}
function saveCurrent(){history.push(makeRecord(exp));renderHistory();drawMorph();}
function restart(){readControls();exp=makeExperiment(config);running=false;$('runBtn').textContent='▶ Run';renderAll();}
function renderAll(){renderSimulation();renderMetrics();renderChart();drawMorph();}
function advance(){if(!exp)restart();if(exp.t>=TMAX){running=false;$('runBtn').textContent='▶ Run';return;}step(exp);renderAll();}
function animation(now){const dt=Math.min(now-lastNow,100);lastNow=now;if(running){acc+=dt;let steps=0;while(acc>=25&&steps<10&&exp.t<TMAX){step(exp);acc-=25;steps++;}renderAll();if(exp.t>=TMAX){running=false;$('runBtn').textContent='▶ Run';}}requestAnimationFrame(animation);}
function sweep(){running=false;$('runBtn').textContent='▶ Run';readControls();const baseline={...config};$('sweepBtn').disabled=true;$('sweepBtn').textContent='Computing…';setTimeout(()=>{let idx=0;for(const p of [180,260,340,420,500])for(const k of [-.15,-.07,.01,.09]){const c={...baseline,pore:p,curvature:k,seed:baseline.seed+idx*101};const e=makeExperiment(c);runToEnd(e,220);history.push(makeRecord(e));idx++;}renderHistory();drawMorph();$('sweepBtn').disabled=false;$('sweepBtn').textContent='Run parameter sweep';$('statusPill').textContent='20 runs saved';},20);}
function csvEscape(v){return '"'+String(v).replaceAll('"','""')+'"';}
function download(filename,body,type){const a=document.createElement('a');const blob=new Blob([body],{type});const url=URL.createObjectURL(blob);a.href=url;a.download=filename;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function csv(){const cols=['id','date','t','architecture','lineage','pore','curvature','stiffness','adhesion','motility','division','chem','seed','population','coverage','commitment','infiltration','matrix','stress','meets_thresholds'];const rows=history.map(r=>[r.id,r.date,r.t,r.c.geometry,r.c.lineage,r.c.pore,r.c.curvature,r.c.stiffness,r.c.adhesion,r.c.motility,r.c.division,r.c.chem,r.c.seed,r.o.population,round(r.o.coverage,3),round(r.o.commitment,3),round(r.o.infiltration,3),round(r.o.matrix,3),round(r.o.stress,3),selectedStatus(r.o)]);download('B4_morphospace_results.csv',cols.join(',')+'\n'+rows.map(r=>r.map(csvEscape).join(',')).join('\n'),'text/csv;charset=utf-8');}
function report(){const records=history.map((r,i)=>`${i+1}. ${r.c.geometry}: pore ${r.c.pore} µm; κ=${r.c.curvature}; E=${r.c.stiffness} GPa; lineage=${r.c.lineage}; cells=${r.o.population}; coverage=${round(r.o.coverage)}%; commitment=${round(r.o.commitment)}%; infiltration=${round(r.o.infiltration)}%; meets criteria=${selectedStatus(r.o)}`).join('\n');const body=`B4 VIRTUAL MORPHOGENESIS LABORATORY — STUDENT REPORT\nGenerated ${new Date().toLocaleString()}\n\nHYPOTHESIS\n${$('hypothesis').value||'(not entered)'}\n\nMETHOD\nOriginal stochastic cell-agent simulator in parameterized scaffold space. No calibration to experimental data. Seeded runs are reproducible within this implementation.\n\nVIABILITY CRITERIA (educational thresholds)\nCoverage >= ${thresholds.coverage}%; committed cells >= ${thresholds.commitment}%; infiltration >= ${thresholds.infiltration}%.\n\nSAVED EXPERIMENTS\n${records||'(none saved)'}\n\nINTERPRETATION AND LIMITATIONS\n${$('notes').value||'(not entered)'}\n\nSCIENTIFIC NOTE\nThe plotted viability heuristic and selected outcome thresholds are pedagogical, not a mathematically established viability kernel or biologically validated design envelope. The model is inspired by, not part of, the published manuscript.\n`;download('B4_experiment_report.txt',body,'text/plain;charset=utf-8');}
const guide=`<p><strong>Purpose.</strong> Investigate how spatial scaffolds and cell-level rules influence emergent organization.</p><ol><li>Select an experiment preset and set pore size, curvature, stiffness, architecture and cell properties.</li><li>Click <strong>Run</strong> or <strong>Step</strong>. Cells migrate, divide, change phenotype and deposit an abstract matrix signal.</li><li>Save each experiment to add it to morphospace. Change one parameter, restart, and save again.</li><li>Use <strong>Run parameter sweep</strong> for twenty simulated designs across pore size and curvature.</li><li>Adjust the three outcome thresholds and ask which saved designs meet them.</li><li>Write your hypothesis and interpretation, then export CSV or a report.</li></ol><p><strong>Recommended activity:</strong> hold everything fixed except pore diameter; compare at least three completed trials. Repeat at a second seed to observe stochastic variation. Ask whether similar coverage can arise from distinct geometric configurations.</p><p><strong>Paper connection:</strong> Figure 6 proposes a two-dimensional scaffold morphospace and Figure 7 illustrates a viability kernel. The simulator offers an independent, simplified teaching analogue, not a reconstruction of those figures' data.</p><div class="figure-links"><a href="assets/paper-figure-6.png" target="_blank" rel="noopener">View manuscript Figure 6</a> · <a href="assets/paper-figure-7.png" target="_blank" rel="noopener">View manuscript Figure 7</a></div>`;
const assumptions=`<p><strong>Scope.</strong> A deliberately reduced, two-dimensional stochastic, off-lattice cell-agent simulation. Scaffold geometry is a periodic obstacle mask. Individual cells experience persistence-free random movement plus directional and neighbour forces, a synthetic geometry-response score, crowding-limited division, a probabilistic phenotype switch and illustrative matrix deposition. These are <strong>hypothetical teaching rules</strong>.</p><p><strong>Not included:</strong> calibrated mechanotransduction, molecular signalling pathways, solvent flow, oxygen diffusion, finite-element stress analysis, realistic pore-scale units, or physical cell deformation. Stiffness modifies synthetic rates and does not solve elasticity equations. Pore diameter and curvature are labels on a parameterized stylized geometry, not a metrologically exact CAD scaffold. The biochemical gradient is a simplified directional cue, not a PDE solution.</p><p><strong>Viability and inference.</strong> The teal ellipse in morphospace is a visual heuristic only. Actual tested points come from simulated trials. Classifications depend on adjustable coverage, differentiation and infiltration thresholds; they are not a mathematical control-theoretic viability kernel. Cells, matrix deposition and fate are synthetic observations, not measurements from the article.</p><p><strong>Reproducibility.</strong> The pseudo-random generator is seeded. Runs are deterministic for fixed parameters and seed within this implementation, but outcomes can change under different seeds and parameter choices.</p><p><strong>Course reference.</strong> Galli, Colangelo, Meleti &amp; Guizzardi, <em>Morphospace Engineering: Morphological Computation in Scaffold Design</em>, <em>BioSystems</em>. The source discusses CA/ABM concepts as conceptual analogies rather than executed simulations; this application is a separate educational implementation.</p>`;
function dialog(title,html){$('dialogTitle').textContent=title;$('dialogText').innerHTML=html;$('infoDialog').showModal();}
function init(){setupControls();$('presetDesc').textContent=presets.colonization.desc;$('preset').onchange=()=>{const p=presets[$('preset').value];config.lineage=p.lineage;config.pore=p.pore;config.curvature=p.curvature;writeControls();$('presetDesc').textContent=p.desc;restart();};
 $('architecture').onchange=()=>restart();$('lineage').onchange=()=>restart();$('chemEnabled').onchange=()=>restart();$('cellCount').onchange=()=>restart();$('seed').onchange=()=>restart();
 $('runBtn').onclick=()=>{if(exp.t>=TMAX)restart();running=!running;$('runBtn').textContent=running?'Ⅱ Pause':'▶ Run';renderMetrics();};$('stepBtn').onclick=()=>{running=false;$('runBtn').textContent='▶ Run';advance();};$('restartBtn').onclick=restart;$('resetBtn').onclick=()=>{Object.assign(config,{pore:300,curvature:-.07,stiffness:1,adhesion:.46,motility:.68,division:.56,geometry:'lattice',lineage:'osteogenic',chem:false,count:65,seed:42});writeControls();restart();};$('saveRunBtn').onclick=saveCurrent;$('sweepBtn').onclick=sweep;$('clearHistoryBtn').onclick=()=>{history=[];renderHistory();drawMorph();};$('csvBtn').onclick=csv;$('reportBtn').onclick=report;$('helpBtn').onclick=()=>dialog('Using the virtual laboratory',guide);$('aboutBtn').onclick=()=>dialog('Model assumptions and limitations',assumptions);$('closeDialog').onclick=()=>$('infoDialog').close();
 $('showPaths').onchange=renderSimulation;$('showField').onchange=renderSimulation;sliders.forEach(([key])=>$('s_'+key).addEventListener('change',restart));
 writeControls();restart();requestAnimationFrame(animation);
}
init();
})();
