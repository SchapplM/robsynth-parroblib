% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR7V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:02:57
% EndTime: 2020-08-07 09:03:00
% DurationCPUTime: 2.83s
% Computational Cost: add. (1530->280), mult. (2595->441), div. (204->14), fcn. (2214->60), ass. (0->216)
t6125 = sin(qJ(1,3));
t6141 = -pkin(5) - pkin(4);
t6090 = t6125 * t6141;
t6134 = cos(qJ(1,3));
t6132 = cos(qJ(3,3));
t6263 = t6132 * pkin(2);
t6085 = pkin(1) + t6263;
t6133 = cos(qJ(2,3));
t6123 = sin(qJ(3,3));
t6124 = sin(qJ(2,3));
t6218 = t6123 * t6124;
t6158 = pkin(2) * t6218 - t6085 * t6133;
t6274 = t6158 * t6134 + t6090;
t6128 = sin(qJ(1,2));
t6091 = t6128 * t6141;
t6137 = cos(qJ(1,2));
t6135 = cos(qJ(3,2));
t6262 = t6135 * pkin(2);
t6087 = pkin(1) + t6262;
t6136 = cos(qJ(2,2));
t6126 = sin(qJ(3,2));
t6127 = sin(qJ(2,2));
t6216 = t6126 * t6127;
t6157 = pkin(2) * t6216 - t6087 * t6136;
t6273 = t6157 * t6137 + t6091;
t6131 = sin(qJ(1,1));
t6092 = t6131 * t6141;
t6140 = cos(qJ(1,1));
t6138 = cos(qJ(3,1));
t6261 = t6138 * pkin(2);
t6089 = pkin(1) + t6261;
t6139 = cos(qJ(2,1));
t6129 = sin(qJ(3,1));
t6130 = sin(qJ(2,1));
t6214 = t6129 * t6130;
t6156 = pkin(2) * t6214 - t6089 * t6139;
t6272 = t6156 * t6140 + t6092;
t6122 = legFrame(1,2);
t6101 = sin(t6122);
t6104 = cos(t6122);
t6077 = t6104 * g(1) - t6101 * g(2);
t6271 = -g(3) * t6131 + t6077 * t6140;
t6121 = legFrame(2,2);
t6100 = sin(t6121);
t6103 = cos(t6121);
t6076 = t6103 * g(1) - t6100 * g(2);
t6270 = -g(3) * t6128 + t6076 * t6137;
t6120 = legFrame(3,2);
t6099 = sin(t6120);
t6102 = cos(t6120);
t6075 = t6102 * g(1) - t6099 * g(2);
t6269 = -g(3) * t6125 + t6075 * t6134;
t6268 = 0.2e1 * pkin(2);
t6267 = 0.2e1 * t6141;
t6260 = -qJ(3,1) + qJ(1,1);
t6259 = qJ(3,1) + qJ(1,1);
t6258 = -qJ(3,2) + qJ(1,2);
t6257 = qJ(3,2) + qJ(1,2);
t6256 = -qJ(3,3) + qJ(1,3);
t6255 = qJ(3,3) + qJ(1,3);
t6254 = qJ(1,1) + 0.2e1 * qJ(2,1);
t6253 = qJ(1,1) - 0.2e1 * qJ(2,1);
t6252 = qJ(1,2) + 0.2e1 * qJ(2,2);
t6251 = qJ(1,2) - 0.2e1 * qJ(2,2);
t6250 = 0.2e1 * qJ(2,3) + qJ(1,3);
t6249 = -0.2e1 * qJ(2,3) + qJ(1,3);
t6072 = t6099 * g(1) + t6102 * g(2);
t6117 = qJ(2,3) + qJ(3,3);
t6093 = sin(t6117);
t6096 = cos(t6117);
t6027 = -t6072 * t6096 + t6093 * t6269;
t6108 = 0.1e1 / t6123;
t6248 = t6027 * t6108;
t6028 = t6072 * t6093 + t6096 * t6269;
t6247 = t6028 * t6108;
t6073 = t6100 * g(1) + t6103 * g(2);
t6118 = qJ(2,2) + qJ(3,2);
t6094 = sin(t6118);
t6097 = cos(t6118);
t6029 = -t6073 * t6097 + t6094 * t6270;
t6109 = 0.1e1 / t6126;
t6246 = t6029 * t6109;
t6030 = t6073 * t6094 + t6097 * t6270;
t6245 = t6030 * t6109;
t6074 = t6101 * g(1) + t6104 * g(2);
t6119 = qJ(2,1) + qJ(3,1);
t6095 = sin(t6119);
t6098 = cos(t6119);
t6031 = -t6074 * t6098 + t6095 * t6271;
t6110 = 0.1e1 / t6129;
t6244 = t6031 * t6110;
t6032 = t6074 * t6095 + t6098 * t6271;
t6243 = t6032 * t6110;
t6242 = (-t6125 * t6158 + t6134 * t6141) * t6108;
t6241 = (-t6128 * t6157 + t6137 * t6141) * t6109;
t6240 = (-t6131 * t6156 + t6140 * t6141) * t6110;
t6078 = 0.1e1 / (t6133 * pkin(1) + pkin(2) * t6096);
t6173 = g(3) * t6134 + t6075 * t6125;
t6239 = t6173 * t6078;
t6079 = 0.1e1 / (t6136 * pkin(1) + pkin(2) * t6097);
t6172 = g(3) * t6137 + t6076 * t6128;
t6238 = t6172 * t6079;
t6080 = 0.1e1 / (t6139 * pkin(1) + pkin(2) * t6098);
t6171 = g(3) * t6140 + t6077 * t6131;
t6237 = t6171 * t6080;
t6236 = 0.1e1 / t6158 * t6108;
t6235 = 0.1e1 / t6157 * t6109;
t6234 = 0.1e1 / t6156 * t6110;
t6230 = t6078 * t6125;
t6229 = t6078 * t6134;
t6228 = t6079 * t6128;
t6227 = t6079 * t6137;
t6226 = t6080 * t6131;
t6225 = t6080 * t6140;
t6111 = t6132 ^ 2;
t6081 = pkin(1) * t6132 + t6111 * t6268 - pkin(2);
t6224 = t6081 * t6134;
t6113 = t6135 ^ 2;
t6082 = pkin(1) * t6135 + t6113 * t6268 - pkin(2);
t6223 = t6082 * t6137;
t6115 = t6138 ^ 2;
t6083 = t6138 * pkin(1) + t6115 * t6268 - pkin(2);
t6222 = t6083 * t6140;
t6221 = (pkin(1) + 0.2e1 * t6263) * t6123;
t6220 = (pkin(1) + 0.2e1 * t6262) * t6126;
t6219 = (pkin(1) + 0.2e1 * t6261) * t6129;
t6217 = t6124 * t6081;
t6215 = t6127 * t6082;
t6213 = t6130 * t6083;
t6212 = t6123 * t6263;
t6211 = t6126 * t6262;
t6210 = t6129 * t6261;
t6209 = t6027 * t6236;
t6208 = t6028 * t6236;
t6207 = t6029 * t6235;
t6206 = t6030 * t6235;
t6205 = t6031 * t6234;
t6204 = t6032 * t6234;
t6033 = -t6072 * t6133 + t6124 * t6269;
t6203 = t6033 * t6236;
t6034 = t6072 * t6124 + t6133 * t6269;
t6202 = t6034 * t6236;
t6035 = -t6073 * t6136 + t6127 * t6270;
t6201 = t6035 * t6235;
t6036 = t6073 * t6127 + t6136 * t6270;
t6200 = t6036 * t6235;
t6037 = -t6074 * t6139 + t6130 * t6271;
t6199 = t6037 * t6234;
t6038 = t6074 * t6130 + t6139 * t6271;
t6198 = t6038 * t6234;
t6197 = t6093 * t6239;
t6196 = t6096 * t6239;
t6195 = t6099 * t6230;
t6194 = t6102 * t6230;
t6193 = t6124 * t6239;
t6192 = t6133 * t6239;
t6191 = t6094 * t6238;
t6190 = t6097 * t6238;
t6189 = t6100 * t6228;
t6188 = t6103 * t6228;
t6187 = t6127 * t6238;
t6186 = t6136 * t6238;
t6185 = t6095 * t6237;
t6184 = t6098 * t6237;
t6183 = t6101 * t6226;
t6182 = t6104 * t6226;
t6181 = t6130 * t6237;
t6180 = t6139 * t6237;
t6179 = t6134 * t6218;
t6178 = t6137 * t6216;
t6177 = t6140 * t6214;
t6142 = 0.2e1 * qJ(3,3);
t6176 = ((cos(qJ(2,3) - t6256) + cos(qJ(2,3) + t6255)) * t6267 + (-sin(0.2e1 * qJ(3,3) - t6249) + sin(t6142 + t6250) + 0.2e1 * t6125) * pkin(2) + (-sin(qJ(3,3) - t6249) + sin(qJ(3,3) + t6250) + sin(t6256) + sin(t6255)) * pkin(1)) / ((-sin(t6142 + qJ(2,3)) + t6124) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t6093) * pkin(1)) / 0.2e1;
t6145 = 0.2e1 * qJ(3,2);
t6175 = ((cos(qJ(2,2) - t6258) + cos(qJ(2,2) + t6257)) * t6267 + (-sin(0.2e1 * qJ(3,2) - t6251) + sin(t6145 + t6252) + 0.2e1 * t6128) * pkin(2) + (-sin(qJ(3,2) - t6251) + sin(qJ(3,2) + t6252) + sin(t6258) + sin(t6257)) * pkin(1)) / ((-sin(t6145 + qJ(2,2)) + t6127) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t6094) * pkin(1)) / 0.2e1;
t6148 = 0.2e1 * qJ(3,1);
t6174 = ((cos(qJ(2,1) - t6260) + cos(qJ(2,1) + t6259)) * t6267 + (-sin(0.2e1 * qJ(3,1) - t6253) + sin(t6148 + t6254) + 0.2e1 * t6131) * pkin(2) + (-sin(qJ(3,1) - t6253) + sin(qJ(3,1) + t6254) + sin(t6260) + sin(t6259)) * pkin(1)) / ((-sin(t6148 + qJ(2,1)) + t6130) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t6095) * pkin(1)) / 0.2e1;
t6170 = t6125 * t6193;
t6169 = t6125 * t6192;
t6168 = t6128 * t6187;
t6167 = t6128 * t6186;
t6166 = t6131 * t6181;
t6165 = t6131 * t6180;
t6164 = t6125 * t6197;
t6163 = t6125 * t6196;
t6162 = t6128 * t6191;
t6161 = t6128 * t6190;
t6160 = t6131 * t6185;
t6159 = t6131 * t6184;
t6155 = pkin(1) * t6179 + (t6179 * t6268 + t6090) * t6132;
t6154 = pkin(1) * t6178 + (t6178 * t6268 + t6091) * t6135;
t6153 = pkin(1) * t6177 + (t6177 * t6268 + t6092) * t6138;
t6152 = 0.1e1 / pkin(1);
t6151 = 0.1e1 / pkin(2);
t6116 = t6139 ^ 2;
t6114 = t6136 ^ 2;
t6112 = t6133 ^ 2;
t6065 = t6129 * pkin(2) * t6139 + t6130 * t6089;
t6064 = t6126 * pkin(2) * t6136 + t6127 * t6087;
t6063 = t6123 * pkin(2) * t6133 + t6124 * t6085;
t6047 = -t6092 * t6214 + (t6115 - 0.1e1) * t6140 * pkin(2);
t6046 = -t6091 * t6216 + (t6113 - 0.1e1) * t6137 * pkin(2);
t6045 = -t6090 * t6218 + (t6111 - 0.1e1) * t6134 * pkin(2);
t6026 = -t6101 * t6065 + t6272 * t6104;
t6025 = -t6104 * t6065 - t6272 * t6101;
t6024 = -t6100 * t6064 + t6273 * t6103;
t6023 = -t6103 * t6064 - t6273 * t6100;
t6022 = -t6099 * t6063 + t6274 * t6102;
t6021 = -t6102 * t6063 - t6274 * t6099;
t6017 = (t6101 * t6219 + t6104 * t6222) * t6116 + (t6101 * t6213 - t6153 * t6104) * t6139 - t6047 * t6104 - t6101 * t6210;
t6016 = (-t6101 * t6222 + t6104 * t6219) * t6116 + (t6153 * t6101 + t6104 * t6213) * t6139 + t6047 * t6101 - t6104 * t6210;
t6015 = (t6100 * t6220 + t6103 * t6223) * t6114 + (t6100 * t6215 - t6154 * t6103) * t6136 - t6046 * t6103 - t6100 * t6211;
t6014 = (-t6100 * t6223 + t6103 * t6220) * t6114 + (t6154 * t6100 + t6103 * t6215) * t6136 + t6046 * t6100 - t6103 * t6211;
t6013 = (t6099 * t6221 + t6102 * t6224) * t6112 + (t6099 * t6217 - t6155 * t6102) * t6133 - t6045 * t6102 - t6099 * t6212;
t6012 = (-t6099 * t6224 + t6102 * t6221) * t6112 + (t6155 * t6099 + t6102 * t6217) * t6133 + t6045 * t6099 - t6102 * t6212;
t1 = [(-t6171 * t6182 - t6172 * t6188 - t6173 * t6194) * MDP(2) + (-t6182 * t6271 - t6188 * t6270 - t6194 * t6269) * MDP(3) + (-t6102 * t6169 - t6103 * t6167 - t6104 * t6165) * MDP(9) + (t6102 * t6170 + t6103 * t6168 + t6104 * t6166) * MDP(10) + (-t6102 * t6163 - t6103 * t6161 - t6104 * t6159) * MDP(16) + (t6102 * t6164 + t6103 * t6162 + t6104 * t6160) * MDP(17) - g(1) * MDP(18) + ((-t6013 * t6203 - t6015 * t6201 - t6017 * t6199) * MDP(9) + (-t6013 * t6202 - t6015 * t6200 - t6017 * t6198) * MDP(10) + (-t6013 * t6209 - t6015 * t6207 - t6017 * t6205) * MDP(16) + (-t6013 * t6208 - t6015 * t6206 - t6017 * t6204) * MDP(17) + ((t6022 * t6248 + t6024 * t6246 + t6026 * t6244) * MDP(16) + (t6022 * t6247 + t6024 * t6245 + t6026 * t6243) * MDP(17)) * t6151) * t6152; (t6171 * t6183 + t6172 * t6189 + t6173 * t6195) * MDP(2) + (t6183 * t6271 + t6189 * t6270 + t6195 * t6269) * MDP(3) + (t6099 * t6169 + t6100 * t6167 + t6101 * t6165) * MDP(9) + (-t6099 * t6170 - t6100 * t6168 - t6101 * t6166) * MDP(10) + (t6099 * t6163 + t6100 * t6161 + t6101 * t6159) * MDP(16) + (-t6099 * t6164 - t6100 * t6162 - t6101 * t6160) * MDP(17) - g(2) * MDP(18) + ((-t6012 * t6203 - t6014 * t6201 - t6016 * t6199) * MDP(9) + (-t6012 * t6202 - t6014 * t6200 - t6016 * t6198) * MDP(10) + (-t6012 * t6209 - t6014 * t6207 - t6016 * t6205) * MDP(16) + (-t6012 * t6208 - t6014 * t6206 - t6016 * t6204) * MDP(17) + ((t6021 * t6248 + t6023 * t6246 + t6025 * t6244) * MDP(16) + (t6021 * t6247 + t6023 * t6245 + t6025 * t6243) * MDP(17)) * t6151) * t6152; (-t6171 * t6225 - t6172 * t6227 - t6173 * t6229) * MDP(2) + (-t6225 * t6271 - t6227 * t6270 - t6229 * t6269) * MDP(3) + (-t6134 * t6192 - t6137 * t6186 - t6140 * t6180 + (t6033 * t6176 + t6035 * t6175 + t6037 * t6174) * t6152) * MDP(9) + (t6134 * t6193 + t6137 * t6187 + t6140 * t6181 + (t6034 * t6176 + t6036 * t6175 + t6038 * t6174) * t6152) * MDP(10) + (-t6134 * t6196 - t6137 * t6190 - t6140 * t6184 + (t6031 * t6174 + t6029 * t6175 + t6027 * t6176 + (t6027 * t6242 + t6029 * t6241 + t6031 * t6240) * t6151) * t6152) * MDP(16) + (t6134 * t6197 + t6137 * t6191 + t6140 * t6185 + (t6032 * t6174 + t6030 * t6175 + t6028 * t6176 + (t6028 * t6242 + t6030 * t6241 + t6032 * t6240) * t6151) * t6152) * MDP(17) - g(3) * MDP(18);];
taugX  = t1;
