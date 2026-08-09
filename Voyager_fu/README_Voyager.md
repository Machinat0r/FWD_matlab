# Voyager 1/2 数据下载程序

本目录提供 MATLAB + Python 下载程序，原始文件来自 NASA/SPDF、CDAWeb
和 JPL/NAIF。数据按卫星、仪器、时间分辨率、质量等级和年份分层保存，
目录形式与现有 MMS 原始数据库相近。默认根目录为：

```text
Z:\SPART-WORK\Data\Voyager
```

## 数据产品

| 程序名 | 内容 | 时间分辨率 | 主要覆盖期 |
|---|---|---:|---|
| `coho1hr` | 合并磁场、等离子体和位置 | 1 h | 任务开始至各产品公开截止时间 |
| `position1day` | HelioWeb 日心位置 | 1 d | 任务开始至 2099 |
| `mag2s` | 早期 MAG 校准数据 | 1.92 s | V1: 1977-09-05 至 1991-12-27；V2: 1977-08-20 至 1991-01-01 |
| `mag2s_unreviewed` | 后期 primary + secondary MAG 实验性数据 | 1.92 s | V1: 1991–1995、2011–2017；V2: 1991–1995、2010–2015、2017–2018 |
| `mag48s_vim` | 后期经审校的 VIM MAG | 48 s | V1: 2009–2025；V2: 2009–2022，目录会实时发现后续更新 |
| `plasma_fine` | PLS Level-3 拟合参数 | 12–192 s | V1: 1977–1980；V2 太阳风: 1977–2007；V2 日球层鞘: 2007–2018 |
| `position1hr` | HelioWeb 日心位置 | 1 h | 任务开始至 2099 |
| `spice_spk` | JPL/NAIF 连续轨道核和飞掠重建解 | 连续多项式 | 任务开始至 2100 |
| `mag48s` | 早期 MAG 48 秒校准数据，可选 | 48 s | 与早期 1.92 秒产品相近 |

`both` 表示两个基础产品，`highres` 表示所有推荐高分辨率产品，
`mag2s_unreviewed` 或 `experimental_mag` 表示后期两个传感器的实验性
1.92 秒产品，`all` 会包括以上全部产品和可选的早期 48 秒 MAG。
高分辨率归档脚本默认先下载
SPDF 产品；JPL 网络可达且安装 MICE 后，可把环境变量
`VOYAGER_INCLUDE_SPICE` 设为 `true`，随后运行同一脚本。

## 最高分辨率和测量精度

- `mag2s` 每个 48 秒源记录包含 25 个 1.92 秒平均值。产品经过校准，
  尚未逐文件科学审校。双磁强计比较给出的验证精度约为
  0.02–0.05 nT，预发射目标精度为 0.1 nT；航天器背景场约为
  0.1–0.2 nT。分析时应保留 `dataConfidence`、`magStatus` 和缺测标记。
- `mag2s_unreviewed` 同样为 1.92 秒平均值，同时保存外置
  primary 与内置 secondary 磁强计。NASA 明确标注该集合通常不具备
  science quality；CDF 中的质量标志约定为 0=good、1=bad。该产品用于
  补充覆盖或诊断时，应同时比较两个传感器，并保留原始质量标志。
- `mag48s_vim` 经过 VIM 团队审校。典型 1σ 为
  `BR/BT/BN` 各约 0.02 nT、`F1` 约 0.03 nT。每个 CDF 内的
  `dF/dBR/dBT/dBN` 给出该文件的实际不确定度，定量分析应优先使用这些值。
- `plasma_fine` 的采样间隔随遥测模式在 12–192 秒之间变化。
  官方说明中的典型 1σ 为速度或径向速度小于 0.5%，密度和热速度
  小于 5%。温度由 `T_eV=0.0052*w^2` 或
  `T_K=60.5*w^2` 计算时，相对误差约为热速度相对误差的两倍。
  Voyager 2 日球层鞘 CDF 还直接提供逐记录的不确定度和 `chi2`。
- `position1hr` 的时间步长为 1 小时。HelioWeb 给出的典型位置精度为
  径向约 1%、角度约 0.1°；小时采样提高时间密度，不改变轨道解本身的
  误差等级。
- `spice_spk` 可在核覆盖范围内按任意时刻计算 J2000 位置和速度。
  加载顺序为 broad trajectory、DE440、2022 x2100、行星飞掠重建解，
  重叠时段会优先采用更精细的 Jupiter、SAT427、URA182 和 NEP097 解。
  x2100 内核在 2030-01-01 给出的估计为：V1 径向 1σ 约
  1.42×10^6 km、经纬角 1σ 约 4.72×10^-4° 和 4.26×10^-4°；
  V2 径向 1σ 约 1.07×10^5 km，三轴位置误差椭球约
  6.92×10^4、2.28×10^5、4.48×10^5 km。轨道可按任意细步长插值，
  远期物理位置误差仍会随外推时间增加。

SPDF 的 `hires1991_2030` 目录含有后期 1.92 秒文件，官方说明指出
这些文件通常不具备科学质量。程序把它们放在带
`experimental_unreviewed_post1991` 标识的独立目录中，且未加入
推荐的 `highres` 组合。官方目录本身是稀疏的：V1 仅有 1991–1995、
2011–2017，V2 仅有 1991–1995、2010–2015、2017–2018；
1996–2008 仍没有这类 1.92 秒文件。后期定量弱磁场分析宜优先采用
经审校的 `mag48s_vim`。

## 存储结构

```text
Voyager/
├── voyager1/
│   ├── coho/1hr/l2/merged_mag_plasma/YYYY/MM/*.cdf
│   ├── mag/1.92s/calibrated_unreviewed/YYYY/MM/*.cdf
│   ├── mag/1.92s/experimental_unreviewed_post1991/
│   │   ├── primary/YYYY/MM/*.cdf
│   │   └── secondary/YYYY/MM/*.cdf
│   ├── mag/48s/reviewed_vim/YYYY/*.cdf
│   ├── pls/hires/l3/solar_wind/YYYY/*.cdf
│   └── ephemeris/
│       ├── 1day/l2/helio_position/YYYY/MM/*.cdf
│       ├── 1hr/l2/helio_position/YYYY/MM/*.cdf
│       └── spice/kernels/...
├── voyager2/
│   ├── coho/1hr/l2/merged_mag_plasma/YYYY/MM/*.cdf
│   ├── mag/1.92s/calibrated_unreviewed/YYYY/MM/*.cdf
│   ├── mag/1.92s/experimental_unreviewed_post1991/
│   │   ├── primary/YYYY/MM/*.cdf
│   │   └── secondary/YYYY/MM/*.cdf
│   ├── mag/48s/reviewed_vim/YYYY/*.cdf
│   ├── pls/hires/l3/{solar_wind,heliosheath}/YYYY/*.cdf
│   └── ephemeris/...
├── voyager_download_manifest.json
├── voyager_highres_download_manifest.json
└── voyager_mag2s_unreviewed_post1991_manifest.json
```

程序先把文件下载到本地暂存目录，核对 HTTP 长度和文件签名后，再通过
目标目录中的临时文件原子发布。SPICE 文件还会逐字节核对官方
SHA-256。重复运行会跳过大小、签名和已知哈希均有效的文件。每次运行
最多同时请求 5 个文件。

若 Windows 上的 Python HTTPS 连接连续超时，可临时设置环境变量
`VOYAGER_FORCE_POWERSHELL_HTTP=true` 并重复运行。下载器会改用
Windows HTTPS 通道；已验证的文件仍会自动跳过。

## MATLAB 使用

下载基础小时/日数据：

```matlab
run_voyager_download
```

补充推荐的高分辨率 SPDF 产品：

```matlab
run_voyager_highres_download
```

下载后期未审校的 primary + secondary 1.92 秒磁场数据：

```matlab
run_voyager_unreviewed_mag_download
```

按时间、卫星和产品选择：

```matlab
report = Voyager_Download( ...
    'Date', '2009-01-01/2010-12-31', ...
    'Spacecraft', 1, ...
    'Products', 'mag48s_vim,position1hr', ...
    'DataRoot', 'Z:\SPART-WORK\Data\Voyager');
```

只查询文件：

```matlab
[names, info] = VoyagerFilenames( ...
    '1977-08-20/2026-07-24', ...
    'Spacecraft', [1 2], ...
    'Products', 'highres');
```

读取 CDF：

```matlab
M = Voyager_Read_CDF(magFile);
P = Voyager_Read_PLS_Hires(plsFile);
plot(M.time, M.B_abs_nT)
```

使用 SPICE 位置需要安装 NASA NAIF MICE：

```matlab
t = datetime(2026, 7, 24, 'TimeZone', 'UTC');
S = Voyager_SPICE_Position(t, 1);
disp(S.position_j2000_km)
```

## Python 后端

后端仅使用 Python 标准库：

```powershell
python download_voyager_files.py `
  --date 1977-08-20/2026-07-24 `
  --spacecraft 1,2 `
  --product mag2s,mag48s_vim,plasma_fine,position1hr `
  --out Z:\SPART-WORK\Data\Voyager `
  --stage C:\Temp\Voyager_staging_highres `
  --manifest-name voyager_highres_download_manifest.json `
  --threads 5
```

后期未审校磁场数据的独立命令：

```powershell
python download_voyager_files.py `
  --date 1991-01-01/2030-12-31 `
  --spacecraft 1,2 `
  --product mag2s_unreviewed `
  --out Z:\SPART-WORK\Data\Voyager `
  --stage C:\Temp\Voyager_staging_mag2s_unreviewed_post1991 `
  --manifest-name voyager_mag2s_unreviewed_post1991_manifest.json `
  --threads 5
```

## 科学使用提示

- Voyager 1 PLS 在 1980 年后没有可用的连续等离子体矩数据。
- Voyager 2 太阳风精细产品止于 2007 年，日球层鞘拟合产品延续至
  2018-11-05。
- 1989 年后磁场接近仪器零点与航天器背景场量级，极弱场研究应把
  VIM 文件中的 `d*` 不确定度和质量说明带入分析。
- `experimental_unreviewed_post1991` 文件应先按 `dataConfidence`
  和相关状态变量筛选，并用 primary/secondary 交叉检查；不要把其
  1.92 秒采样间隔当成更高的绝对测量精度。
- 整时、整日或整年 CDF 会保留官方补齐记录。两个读取器会依据
  `FILLVAL`、`VALIDMIN` 和 `VALIDMAX` 将无效数据转为 `NaN`。
- 原始 CDF、全局属性、变量属性、版本号和文件名均保留。

## 官方来源

- NASA/SPDF Voyager 1: https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager1/
- NASA/SPDF Voyager 2: https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/
- Voyager 1 MAG 官方质量说明: https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager1/magnetic_fields_cdaweb/00readme.txt
- Voyager 2 MAG 官方质量说明: https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/magnetic_fields_cdaweb/00readme.txt
- CDAWeb Voyager Notes: https://cdaweb.gsfc.nasa.gov/misc/NotesV.html
- JPL/NAIF Voyager kernels: https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/
- Voyager 1 1.92 s MAG DOI: https://doi.org/10.48322/qhrq-xt70
- Voyager 2 1.92 s MAG DOI: https://doi.org/10.48322/d4s2-1q92
- Voyager 1 VIM 48 s MAG DOI: https://doi.org/10.48322/envb-3w78
- Voyager 2 VIM 48 s MAG DOI: https://doi.org/10.48322/kgbb-j403
- Voyager 1 PLS high-resolution DOI: https://doi.org/10.48322/kcn9-y866
- Voyager 2 PLS high-resolution DOI: https://doi.org/10.48322/3r7q-jg89
