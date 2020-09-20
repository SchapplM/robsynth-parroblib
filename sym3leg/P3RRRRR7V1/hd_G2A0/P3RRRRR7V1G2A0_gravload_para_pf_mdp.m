% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V1G2A0
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
%   see P3RRRRR7V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:52:33
% EndTime: 2020-08-07 03:52:36
% DurationCPUTime: 2.89s
% Computational Cost: add. (1530->280), mult. (2595->435), div. (204->14), fcn. (2214->60), ass. (0->216)
t6168 = legFrame(2,2);
t6145 = sin(t6168);
t6148 = cos(t6168);
t6124 = t6148 * g(1) - t6145 * g(2);
t6175 = sin(qJ(1,2));
t6184 = cos(qJ(1,2));
t6315 = -g(3) * t6175 + t6124 * t6184;
t6314 = 2 * pkin(2);
t6188 = (pkin(4) + pkin(5));
t6313 = 2 * t6188;
t6179 = cos(qJ(3,3));
t6311 = t6179 * pkin(2);
t6182 = cos(qJ(3,2));
t6310 = t6182 * pkin(2);
t6185 = cos(qJ(3,1));
t6309 = t6185 * pkin(2);
t6308 = -qJ(3,1) + qJ(1,1);
t6307 = qJ(3,1) + qJ(1,1);
t6306 = -qJ(3,2) + qJ(1,2);
t6305 = qJ(3,2) + qJ(1,2);
t6304 = -qJ(3,3) + qJ(1,3);
t6303 = qJ(3,3) + qJ(1,3);
t6302 = (2 * qJ(2,1)) + qJ(1,1);
t6301 = -(2 * qJ(2,1)) + qJ(1,1);
t6300 = (2 * qJ(2,2)) + qJ(1,2);
t6299 = -(2 * qJ(2,2)) + qJ(1,2);
t6298 = (2 * qJ(2,3)) + qJ(1,3);
t6297 = -(2 * qJ(2,3)) + qJ(1,3);
t6167 = legFrame(3,2);
t6144 = sin(t6167);
t6147 = cos(t6167);
t6120 = t6144 * g(1) + t6147 * g(2);
t6164 = qJ(2,3) + qJ(3,3);
t6138 = sin(t6164);
t6141 = cos(t6164);
t6123 = t6147 * g(1) - t6144 * g(2);
t6172 = sin(qJ(1,3));
t6181 = cos(qJ(1,3));
t6222 = g(3) * t6181 + t6123 * t6172;
t6075 = -t6120 * t6141 + t6138 * t6222;
t6170 = sin(qJ(3,3));
t6155 = 0.1e1 / t6170;
t6296 = t6075 * t6155;
t6076 = t6120 * t6138 + t6141 * t6222;
t6295 = t6076 * t6155;
t6169 = legFrame(1,2);
t6146 = sin(t6169);
t6149 = cos(t6169);
t6122 = t6146 * g(1) + t6149 * g(2);
t6166 = qJ(2,1) + qJ(3,1);
t6140 = sin(t6166);
t6143 = cos(t6166);
t6125 = t6149 * g(1) - t6146 * g(2);
t6178 = sin(qJ(1,1));
t6187 = cos(qJ(1,1));
t6220 = g(3) * t6187 + t6125 * t6178;
t6077 = -t6122 * t6143 + t6140 * t6220;
t6176 = sin(qJ(3,1));
t6157 = 0.1e1 / t6176;
t6294 = t6077 * t6157;
t6078 = t6122 * t6140 + t6143 * t6220;
t6293 = t6078 * t6157;
t6121 = t6145 * g(1) + t6148 * g(2);
t6165 = qJ(2,2) + qJ(3,2);
t6139 = sin(t6165);
t6142 = cos(t6165);
t6221 = g(3) * t6184 + t6124 * t6175;
t6079 = -t6121 * t6142 + t6139 * t6221;
t6173 = sin(qJ(3,2));
t6156 = 0.1e1 / t6173;
t6292 = t6079 * t6156;
t6080 = t6121 * t6139 + t6142 * t6221;
t6291 = t6080 * t6156;
t6133 = pkin(1) + t6311;
t6180 = cos(qJ(2,3));
t6171 = sin(qJ(2,3));
t6268 = t6170 * t6171;
t6205 = pkin(2) * t6268 - t6133 * t6180;
t6290 = (-t6172 * t6188 + t6205 * t6181) * t6155;
t6135 = pkin(1) + t6310;
t6183 = cos(qJ(2,2));
t6174 = sin(qJ(2,2));
t6266 = t6173 * t6174;
t6204 = pkin(2) * t6266 - t6135 * t6183;
t6289 = (-t6175 * t6188 + t6204 * t6184) * t6156;
t6137 = pkin(1) + t6309;
t6186 = cos(qJ(2,1));
t6177 = sin(qJ(2,1));
t6264 = t6176 * t6177;
t6203 = pkin(2) * t6264 - t6137 * t6186;
t6288 = (-t6178 * t6188 + t6203 * t6187) * t6157;
t6105 = -g(3) * t6172 + t6123 * t6181;
t6126 = 0.1e1 / (t6180 * pkin(1) + pkin(2) * t6141);
t6287 = t6105 * t6126;
t6127 = 0.1e1 / (t6183 * pkin(1) + pkin(2) * t6142);
t6286 = t6315 * t6127;
t6106 = -g(3) * t6178 + t6125 * t6187;
t6128 = 0.1e1 / (t6186 * pkin(1) + pkin(2) * t6143);
t6285 = t6106 * t6128;
t6284 = 0.1e1 / t6205 * t6155;
t6283 = 0.1e1 / t6204 * t6156;
t6282 = 0.1e1 / t6203 * t6157;
t6280 = t6126 * t6172;
t6279 = t6126 * t6181;
t6278 = t6127 * t6175;
t6277 = t6127 * t6184;
t6276 = t6128 * t6178;
t6275 = t6128 * t6187;
t6158 = t6179 ^ 2;
t6129 = pkin(1) * t6179 + t6158 * t6314 - pkin(2);
t6274 = t6129 * t6172;
t6160 = t6182 ^ 2;
t6130 = pkin(1) * t6182 + t6160 * t6314 - pkin(2);
t6273 = t6130 * t6175;
t6162 = t6185 ^ 2;
t6131 = t6185 * pkin(1) + t6162 * t6314 - pkin(2);
t6272 = t6131 * t6178;
t6271 = (pkin(1) + 0.2e1 * t6311) * t6170;
t6270 = (pkin(1) + 0.2e1 * t6310) * t6173;
t6269 = (pkin(1) + 0.2e1 * t6309) * t6176;
t6267 = t6171 * t6129;
t6265 = t6174 * t6130;
t6263 = t6177 * t6131;
t6262 = t6181 * t6188;
t6261 = t6184 * t6188;
t6260 = t6187 * t6188;
t6259 = t6170 * t6311;
t6258 = t6173 * t6310;
t6257 = t6176 * t6309;
t6256 = t6075 * t6284;
t6255 = t6076 * t6284;
t6254 = t6077 * t6282;
t6253 = t6078 * t6282;
t6252 = t6079 * t6283;
t6251 = t6080 * t6283;
t6081 = -t6120 * t6180 + t6171 * t6222;
t6250 = t6081 * t6284;
t6082 = t6120 * t6171 + t6180 * t6222;
t6249 = t6082 * t6284;
t6083 = -t6122 * t6186 + t6177 * t6220;
t6248 = t6083 * t6282;
t6084 = t6122 * t6177 + t6186 * t6220;
t6247 = t6084 * t6282;
t6085 = -t6121 * t6183 + t6174 * t6221;
t6246 = t6085 * t6283;
t6086 = t6121 * t6174 + t6183 * t6221;
t6245 = t6086 * t6283;
t6244 = t6141 * t6287;
t6243 = t6144 * t6279;
t6242 = t6147 * t6279;
t6241 = t6180 * t6287;
t6240 = t6139 * t6286;
t6239 = t6142 * t6286;
t6238 = t6145 * t6277;
t6237 = t6148 * t6277;
t6236 = t6174 * t6286;
t6235 = t6183 * t6286;
t6234 = t6143 * t6285;
t6233 = t6146 * t6275;
t6232 = t6149 * t6275;
t6231 = t6186 * t6285;
t6230 = t6172 * t6268;
t6229 = t6105 * t6280;
t6228 = t6175 * t6266;
t6227 = t6178 * t6264;
t6226 = t6106 * t6276;
t6189 = 0.2e1 * qJ(3,3);
t6225 = ((sin(qJ(2,3) - t6304) - sin(qJ(2,3) + t6303)) * t6313 + (-cos(0.2e1 * qJ(3,3) - t6297) - cos(t6189 + t6298) - 0.2e1 * t6181) * pkin(2) + (-cos(qJ(3,3) - t6297) - cos(qJ(3,3) + t6298) - cos(t6304) - cos(t6303)) * pkin(1)) / ((-sin(t6189 + qJ(2,3)) + t6171) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t6138) * pkin(1)) / 0.2e1;
t6192 = 0.2e1 * qJ(3,2);
t6224 = ((sin(qJ(2,2) - t6306) - sin(qJ(2,2) + t6305)) * t6313 + (-cos(0.2e1 * qJ(3,2) - t6299) - cos(t6192 + t6300) - 0.2e1 * t6184) * pkin(2) + (-cos(qJ(3,2) - t6299) - cos(qJ(3,2) + t6300) - cos(t6306) - cos(t6305)) * pkin(1)) / ((-sin(t6192 + qJ(2,2)) + t6174) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t6139) * pkin(1)) / 0.2e1;
t6195 = 0.2e1 * qJ(3,1);
t6223 = ((sin(qJ(2,1) - t6308) - sin(qJ(2,1) + t6307)) * t6313 + (-cos(0.2e1 * qJ(3,1) - t6301) - cos(t6195 + t6302) - 0.2e1 * t6187) * pkin(2) + (-cos(qJ(3,1) - t6301) - cos(qJ(3,1) + t6302) - cos(t6308) - cos(t6307)) * pkin(1)) / ((-sin(t6195 + qJ(2,1)) + t6177) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t6140) * pkin(1)) / 0.2e1;
t6219 = t6181 * t6241;
t6218 = t6184 * t6236;
t6217 = t6184 * t6235;
t6216 = t6187 * t6231;
t6215 = t6105 * t6243;
t6214 = t6106 * t6233;
t6213 = t6105 * t6242;
t6212 = t6106 * t6232;
t6211 = t6105 * t6138 * t6279;
t6210 = t6181 * t6244;
t6209 = t6184 * t6240;
t6208 = t6184 * t6239;
t6207 = t6106 * t6140 * t6275;
t6206 = t6187 * t6234;
t6202 = pkin(1) * t6230 + (t6230 * t6314 + t6262) * t6179;
t6201 = pkin(1) * t6228 + (t6228 * t6314 + t6261) * t6182;
t6200 = pkin(1) * t6227 + (t6227 * t6314 + t6260) * t6185;
t6199 = 0.1e1 / pkin(1);
t6198 = 1 / pkin(2);
t6163 = t6186 ^ 2;
t6161 = t6183 ^ 2;
t6159 = t6180 ^ 2;
t6116 = t6176 * pkin(2) * t6186 + t6177 * t6137;
t6115 = t6173 * pkin(2) * t6183 + t6174 * t6135;
t6114 = t6170 * pkin(2) * t6180 + t6171 * t6133;
t6098 = -t6260 * t6264 + (t6162 - 0.1e1) * t6178 * pkin(2);
t6097 = -t6261 * t6266 + (t6160 - 0.1e1) * t6175 * pkin(2);
t6096 = -t6262 * t6268 + (t6158 - 0.1e1) * t6172 * pkin(2);
t6089 = -t6203 * t6178 - t6260;
t6088 = -t6204 * t6175 - t6261;
t6087 = -t6205 * t6172 - t6262;
t6074 = -t6089 * t6149 - t6116 * t6146;
t6073 = t6089 * t6146 - t6116 * t6149;
t6072 = -t6088 * t6148 - t6115 * t6145;
t6071 = t6088 * t6145 - t6115 * t6148;
t6070 = -t6087 * t6147 - t6114 * t6144;
t6069 = t6087 * t6144 - t6114 * t6147;
t6065 = (t6146 * t6269 + t6149 * t6272) * t6163 + (t6146 * t6263 - t6200 * t6149) * t6186 - t6098 * t6149 - t6146 * t6257;
t6064 = (-t6146 * t6272 + t6149 * t6269) * t6163 + (t6200 * t6146 + t6149 * t6263) * t6186 + t6098 * t6146 - t6149 * t6257;
t6063 = (t6145 * t6270 + t6148 * t6273) * t6161 + (t6145 * t6265 - t6201 * t6148) * t6183 - t6097 * t6148 - t6145 * t6258;
t6062 = (-t6145 * t6273 + t6148 * t6270) * t6161 + (t6201 * t6145 + t6148 * t6265) * t6183 + t6097 * t6145 - t6148 * t6258;
t6061 = (t6144 * t6271 + t6147 * t6274) * t6159 + (t6144 * t6267 - t6202 * t6147) * t6180 - t6096 * t6147 - t6144 * t6259;
t6060 = (-t6144 * t6274 + t6147 * t6271) * t6159 + (t6202 * t6144 + t6147 * t6267) * t6180 + t6096 * t6144 - t6147 * t6259;
t1 = [(-t6237 * t6315 - t6212 - t6213) * MDP(2) + (t6220 * t6232 + t6221 * t6237 + t6222 * t6242) * MDP(3) + (-t6147 * t6219 - t6148 * t6217 - t6149 * t6216) * MDP(9) + (t6148 * t6218 + t6171 * t6213 + t6177 * t6212) * MDP(10) + (-t6147 * t6210 - t6148 * t6208 - t6149 * t6206) * MDP(16) + (t6147 * t6211 + t6148 * t6209 + t6149 * t6207) * MDP(17) - g(1) * MDP(18) + ((-t6061 * t6250 - t6063 * t6246 - t6065 * t6248) * MDP(9) + (-t6061 * t6249 - t6063 * t6245 - t6065 * t6247) * MDP(10) + (-t6061 * t6256 - t6063 * t6252 - t6065 * t6254) * MDP(16) + (-t6061 * t6255 - t6063 * t6251 - t6065 * t6253) * MDP(17) + ((t6070 * t6296 + t6072 * t6292 + t6074 * t6294) * MDP(16) + (t6070 * t6295 + t6072 * t6291 + t6074 * t6293) * MDP(17)) * t6198) * t6199; (t6238 * t6315 + t6214 + t6215) * MDP(2) + (-t6220 * t6233 - t6221 * t6238 - t6222 * t6243) * MDP(3) + (t6144 * t6219 + t6145 * t6217 + t6146 * t6216) * MDP(9) + (-t6145 * t6218 - t6171 * t6215 - t6177 * t6214) * MDP(10) + (t6144 * t6210 + t6145 * t6208 + t6146 * t6206) * MDP(16) + (-t6144 * t6211 - t6145 * t6209 - t6146 * t6207) * MDP(17) - g(2) * MDP(18) + ((-t6060 * t6250 - t6062 * t6246 - t6064 * t6248) * MDP(9) + (-t6060 * t6249 - t6062 * t6245 - t6064 * t6247) * MDP(10) + (-t6060 * t6256 - t6062 * t6252 - t6064 * t6254) * MDP(16) + (-t6060 * t6255 - t6062 * t6251 - t6064 * t6253) * MDP(17) + ((t6069 * t6296 + t6071 * t6292 + t6073 * t6294) * MDP(16) + (t6069 * t6295 + t6071 * t6291 + t6073 * t6293) * MDP(17)) * t6198) * t6199; (t6278 * t6315 + t6226 + t6229) * MDP(2) + (-t6220 * t6276 - t6221 * t6278 - t6222 * t6280) * MDP(3) + (t6172 * t6241 + t6175 * t6235 + t6178 * t6231 + (t6081 * t6225 + t6083 * t6223 + t6085 * t6224) * t6199) * MDP(9) + (-t6175 * t6236 - t6171 * t6229 - t6177 * t6226 + (t6082 * t6225 + t6084 * t6223 + t6086 * t6224) * t6199) * MDP(10) + (t6172 * t6244 + t6175 * t6239 + t6178 * t6234 + (t6077 * t6223 + t6079 * t6224 + t6075 * t6225 + (t6075 * t6290 + t6077 * t6288 + t6079 * t6289) * t6198) * t6199) * MDP(16) + (-t6175 * t6240 - t6138 * t6229 - t6140 * t6226 + (t6078 * t6223 + t6080 * t6224 + t6076 * t6225 + (t6076 * t6290 + t6078 * t6288 + t6080 * t6289) * t6198) * t6199) * MDP(17) - g(3) * MDP(18);];
taugX  = t1;
