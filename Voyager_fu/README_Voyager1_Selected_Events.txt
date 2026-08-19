Voyager 1 圈选事件绘图与数据库
================================

运行入口：
  run_voyager1_selected_events.m

主程序：
  Voyager_Interstellar_Monthly\Voyager1_Plot_Selected_Events.m

默认输入：
  Z:\SPART-WORK\Data\Voyager\voyager1\coho\1hr\l2\merged_mag_plasma

默认图片输出：
  Voyager_Interstellar_Monthly\Voyager1_Selected_Events_1h

默认 SQLite 数据库：
  Z:\SPART-WORK\Data\Voyager\voyager1\events\Voyager1_Selected_Events_1h.sqlite

数据库表：
  events                     PPT 圈选区间、来源页码、原始映射时间、绘图时间
  event_sources              每个事件使用的原始 CDF 文件
  magnetic_field             每个原始时刻的 |B|、BR、BT、BN、日心距离
  particle_flux              LECP 与 CRS 全部能道的长表数据
  provenance                 精度、缺测处理、圆圈坐标映射方法
  literature_style_references 绘图样式参考文章

处理约束：
  1. 只读取仪器已有记录。
  2. 不平滑、不插值、不重采样、不补点。
  3. 缺测在 SQLite 中保存为 NULL；图中断线或留白。
  4. LECP 画三条对数通量曲线；CRS 画时间-能量-对数通量谱图。
  5. COHO CDF 的 CRS 元数据标注为 6 小时产品；程序保留该原始节奏。
