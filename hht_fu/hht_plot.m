function fn = hht_plot(T, F, P, varargin)
%------written by Wending Fu, Nov.2023 in Beijing------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       南无电子阿弥陀佛驱散仿生bug
%                                _ooOoo_
%                               o8888888o
%                               88" . "88
%                               (| -_- |)
%                               O\  =  /O
%                            ____/`---'\____
%                          .'  \\|     |//  `.
%                         /  \\|||  :  |||//  \
%                        /  _||||| -:- |||||-  \
%                        |   | \\\  -  /// |   |
%                        | \_|  ''\-/''  |   |
%                        \  .-\__  `-`  ___/-. /
%                      ___`. .'  /-.-\  `. . __
%                   ."" '<  `.___\_<|>_/___.'  >'"".
%                  | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%                  \  \ `-.   \_ __\ /__ _/   .-` /  /
%             ======`-.____`-.___\_____/___.-`____.-'======
% 	                   `=-='
%                 天地玄宗，万气本根。广修亿劫，证吾神通。
%                 三界内外，惟道独尊。体有金光，覆映吾身。
%                 视之不见，听之不闻。包罗天地，养育群生。
%                 受持万遍，身有光明。三界侍卫，五帝司迎。
%                 万神朝礼，役使雷霆。鬼妖丧胆，精怪忘形。
%                 内有霹雳，雷神隐名。洞慧交彻，五炁腾腾。
%                金光速现，覆护真人。急急如律令，bug全去除！
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input
parser = inputParser;
fs = 1/median(diff(T));
defaultFRange = [1, fs/2];
nanP = P; nanP(P == 0) = nan;
defaultPRange = log10([nanmin(nanP,[],'all'), nanmax(nanP,[],'all')]);
defaultLineWidth = 1;
defaultTimeAxis = 1;
defaultFaceColor = 'none';
defaultEdgeColor = 'interp';
defaultEdgeAlpha = 'interp';
addParameter(parser, 'FRange', defaultFRange, @ismatrix);
addParameter(parser, 'PRange', defaultPRange, @ismatrix);
addParameter(parser, 'LineWidth', defaultLineWidth, @isnumeric);
addParameter(parser, 'TimeAxis', defaultTimeAxis, @isnumeric);
addParameter(parser, 'FaceColor', defaultFaceColor, @ischar);
addParameter(parser, 'EdgeColor', defaultEdgeColor, @ischar);
addParameter(parser, 'EdgeAlpha', defaultEdgeAlpha);

parse(parser, varargin{:});

FRange = parser.Results.FRange;
PRange = parser.Results.PRange;
LineWidth = parser.Results.LineWidth;
TimeAxisFlag = parser.Results.TimeAxis;
FaceColor = parser.Results.FaceColor;
EdgeColor = parser.Results.EdgeColor;
EdgeAlpha = parser.Results.EdgeAlpha;
switch TimeAxisFlag
case 0
    TimeAxis = 'nolabels';
case 1
    TimeAxis = 'date';
case 2
    TimeAxis = 'nodate';
end
%% plot
if size(F,2) ~= size(P,2), warning('wrong dimension of F and P'); help hht_plot; end
for imf = 1:size(F,2)
tempF = F(:,imf);
tempP = P(:,imf);
tempP_log10 = log10(tempP);
fn = patch([T(1);T;T(end)], [0;tempF;0], [nan;tempP_log10;nan], ...
    'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,...
    'FaceColor', FaceColor, 'FaceVertexAlphaData',[nan;tempP;nan],...
    'LineWidth', LineWidth, 'FaceAlpha', 'interp');
hold on;
end

xyrange = [0,T(end),FRange(1),FRange(2)];
axis(xyrange);
set(gca, 'yscale','log');
irf_timeaxis(gca,TimeAxis);
clim(PRange);
colormap jet
colorbar