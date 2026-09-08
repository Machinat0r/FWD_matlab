import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
const out='Z:/SPART-WORK/Data/MMS/derived/KH/catalog_audit_20260903';
const temp=path.join(os.tmpdir(),'KH_catalog_20260903');
const parse=s=>JSON.parse(s.replace(/^\uFEFF/,''));
const input=parse(await fs.readFile(path.join(temp,'input_bundle.json'),'utf8'));
const events=[...input.catalog,{EventID:'KH090',StartUTC:'2018-11-06 13:00:00',EndUTC:'2018-11-06 14:10:00'},
{EventID:'CAND01',StartUTC:'2015-10-02 08:20:00',EndUTC:'2015-10-02 11:00:00'}];
const specs=[
['B','brst','fgm',''],['Vi','brst','fpi','dis-moms'],['Ve','brst','fpi','des-moms'],['E','brst','edp','dce'],
['B','srvy','fgm',''],['Vi','fast','fpi','dis-moms'],['Ve','fast','fpi','des-moms'],['E','fast','edp','dce']];
const cacheDir=path.join(out,'sdc_metadata');await fs.mkdir(cacheDir,{recursive:true});
const days=e=>{const a=[];for(let d=new Date(e.StartUTC.slice(0,10));d<=new Date(e.EndUTC.slice(0,10));d.setUTCDate(d.getUTCDate()+1))a.push(d.toISOString().slice(0,10));return a};
const keys=[...new Set(events.flatMap(e=>days(e)))];
const requests=keys.flatMap(day=>specs.slice(0,4).map(s=>({day,s})));
const metadata=new Map();
async function get(day,s){
 const key=[day,...s].join('_');if(metadata.has(key))return metadata.get(key);
 const filename=path.join(cacheDir,key+'.json');
 let data;try{data=parse(await fs.readFile(filename,'utf8'))}catch{}
 if(!data){data={url:'',error:'缺少SDC缓存',files:[]}}
 if(false){
  const end=new Date(day);end.setUTCDate(end.getUTCDate()+1);
  const u=new URL('https://lasp.colorado.edu/mms/sdc/public/files/api/v1/file_info/science');
  const params={sc_id:'mms1,mms2,mms3,mms4',instrument_id:s[2],data_rate_mode:s[1],data_level:'l2',start_date:day,end_date:end.toISOString().slice(0,10)};
  if(s[3])params.descriptor=s[3];u.search=new URLSearchParams(params).toString();
  data={url:u.href,checkedUTC:new Date().toISOString(),error:'',files:[]};
  for(let attempt=0;attempt<2;attempt++){
   try{const r=await fetch(u,{signal:AbortSignal.timeout(30000)});if(!r.ok)throw Error('HTTP '+r.status);
    const x=await r.json();data={...data,...x,error:''};break;
   }catch(e){data.error=String(e)}
  }
  await fs.writeFile(filename,JSON.stringify(data));
 }
 metadata.set(key,data);return data;
}
let next=0,done=0;
await Promise.all(Array.from({length:4},async()=>{while(next<requests.length){const {day,s}=requests[next++];await get(day,s);done++;if(done%40===0)console.log('SDC burst requests '+done+'/'+requests.length)}}));
async function summarize(e,s,scope){
 const ds=await Promise.all(days(e).map(d=>get(d,s)));
 const files=ds.flatMap(x=>x.files||[]);
 const t1=Date.parse(e.StartUTC.replace(' ','T')+'Z'),t2=Date.parse(e.EndUTC.replace(' ','T')+'Z');
 const selected=files.filter(f=>{if(scope==='当日文件')return days(e).includes(f.timetag.slice(0,10));const t=Date.parse(f.timetag+'Z');return t>=t1&&t<t2});
 const counts=[1,2,3,4].map(ic=>selected.filter(f=>f.file_name.startsWith('mms'+ic+'_')).length);
 return {EventID:e.EventID,Product:s[0],Mode:s[1],Scope:scope,Counts:counts,SpacecraftCount:counts.filter(n=>n>0).length,
  Files:selected,URLs:ds.map(x=>x.url),Truncated:ds.some(x=>x.files_were_truncated),Errors:ds.map(x=>x.error).filter(Boolean)};
}
const rows=[];
for(const e of events){
 const burst=[];for(const s of specs.slice(0,4)){const r=await summarize(e,s,'文件起始时刻在事件内');rows.push(r);burst.push(r)}
 if(burst.some(r=>r.SpacecraftCount<4||r.Truncated||r.Errors.length)||e.EventID.startsWith('CAND'))
  for(const s of specs.slice(4))rows.push(await summarize(e,s,'当日文件'));
}
await fs.writeFile(path.join(out,'public_availability.json'),JSON.stringify({checkedUTC:new Date().toISOString(),events:events.length,rows,complete:true}));
console.log('Finished public metadata for '+events.length+' events; '+rows.length+' product rows.');
