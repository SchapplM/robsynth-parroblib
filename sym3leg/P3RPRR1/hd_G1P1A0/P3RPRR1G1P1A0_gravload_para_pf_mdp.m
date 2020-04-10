% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:23
% EndTime: 2020-03-09 21:23:24
% DurationCPUTime: 0.42s
% Computational Cost: add. (535->86), mult. (435->127), div. (54->4), fcn. (402->42), ass. (0->78)
t1313 = MDP(4) * pkin(1) + MDP(2);
t1282 = legFrame(3,3);
t1311 = cos(t1282);
t1283 = legFrame(2,3);
t1310 = cos(t1283);
t1284 = legFrame(1,3);
t1309 = cos(t1284);
t1273 = t1282 + qJ(1,3);
t1274 = t1283 + qJ(1,2);
t1275 = t1284 + qJ(1,1);
t1294 = pkin(7) + qJ(3,3);
t1249 = 0.1e1 / (pkin(1) * sin(t1294) + sin(qJ(3,3)) * pkin(2));
t1261 = pkin(7) + t1273;
t1258 = qJ(3,3) + t1261;
t1252 = sin(t1258);
t1308 = (-pkin(1) * sin(t1273) - pkin(2) * sin(t1261) - pkin(3) * t1252) * t1249;
t1295 = pkin(7) + qJ(3,2);
t1250 = 0.1e1 / (pkin(1) * sin(t1295) + sin(qJ(3,2)) * pkin(2));
t1262 = pkin(7) + t1274;
t1259 = qJ(3,2) + t1262;
t1253 = sin(t1259);
t1307 = (-pkin(1) * sin(t1274) - pkin(2) * sin(t1262) - pkin(3) * t1253) * t1250;
t1296 = pkin(7) + qJ(3,1);
t1251 = 0.1e1 / (pkin(1) * sin(t1296) + sin(qJ(3,1)) * pkin(2));
t1263 = pkin(7) + t1275;
t1260 = qJ(3,1) + t1263;
t1254 = sin(t1260);
t1306 = (-pkin(1) * sin(t1275) - pkin(2) * sin(t1263) - pkin(3) * t1254) * t1251;
t1255 = cos(t1258);
t1305 = (-pkin(1) * cos(t1273) - pkin(2) * cos(t1261) - pkin(3) * t1255) * t1249;
t1256 = cos(t1259);
t1304 = (-pkin(1) * cos(t1274) - pkin(2) * cos(t1262) - pkin(3) * t1256) * t1250;
t1257 = cos(t1260);
t1303 = (-pkin(1) * cos(t1275) - pkin(2) * cos(t1263) - pkin(3) * t1257) * t1251;
t1302 = t1249 * t1252;
t1301 = t1249 * t1255;
t1300 = t1250 * t1253;
t1299 = t1250 * t1256;
t1298 = t1251 * t1254;
t1297 = t1251 * t1257;
t1291 = 0.1e1 / pkin(3);
t1290 = cos(qJ(1,1));
t1289 = cos(qJ(1,2));
t1288 = cos(qJ(1,3));
t1287 = sin(qJ(1,1));
t1286 = sin(qJ(1,2));
t1285 = sin(qJ(1,3));
t1281 = sin(t1284);
t1280 = sin(t1283);
t1279 = sin(t1282);
t1278 = qJ(1,1) + t1296;
t1277 = qJ(1,2) + t1295;
t1276 = qJ(1,3) + t1294;
t1269 = cos(t1278);
t1268 = cos(t1277);
t1267 = cos(t1276);
t1266 = sin(t1278);
t1265 = sin(t1277);
t1264 = sin(t1276);
t1248 = g(1) * t1309 + t1281 * g(2);
t1247 = g(1) * t1310 + t1280 * g(2);
t1246 = g(1) * t1311 + t1279 * g(2);
t1245 = t1281 * g(1) - g(2) * t1309;
t1244 = t1280 * g(1) - g(2) * t1310;
t1243 = t1279 * g(1) - g(2) * t1311;
t1236 = -t1245 * t1287 + t1248 * t1290;
t1235 = -t1244 * t1286 + t1247 * t1289;
t1234 = -t1243 * t1285 + t1246 * t1288;
t1233 = t1245 * t1290 + t1248 * t1287;
t1232 = t1244 * t1289 + t1247 * t1286;
t1231 = t1243 * t1288 + t1246 * t1285;
t1230 = -t1245 * t1266 + t1248 * t1269;
t1229 = -t1244 * t1265 + t1247 * t1268;
t1228 = -t1243 * t1264 + t1246 * t1267;
t1227 = t1245 * t1269 + t1248 * t1266;
t1226 = t1244 * t1268 + t1247 * t1265;
t1225 = t1243 * t1267 + t1246 * t1264;
t1 = [(t1234 * t1301 + t1235 * t1299 + t1236 * t1297) * MDP(3) + (t1225 * t1301 + t1226 * t1299 + t1227 * t1297) * MDP(6) + (t1228 * t1301 + t1229 * t1299 + t1230 * t1297) * MDP(7) - g(1) * MDP(8) + ((t1225 * t1305 + t1226 * t1304 + t1227 * t1303) * MDP(6) + (t1228 * t1305 + t1229 * t1304 + t1230 * t1303) * MDP(7)) * t1291 + t1313 * (t1231 * t1301 + t1232 * t1299 + t1233 * t1297); (t1234 * t1302 + t1235 * t1300 + t1236 * t1298) * MDP(3) + (t1225 * t1302 + t1226 * t1300 + t1227 * t1298) * MDP(6) + (t1228 * t1302 + t1229 * t1300 + t1230 * t1298) * MDP(7) - g(2) * MDP(8) + ((t1225 * t1308 + t1226 * t1307 + t1227 * t1306) * MDP(6) + (t1228 * t1308 + t1229 * t1307 + t1230 * t1306) * MDP(7)) * t1291 + t1313 * (t1231 * t1302 + t1232 * t1300 + t1233 * t1298); (-(3 * MDP(4)) - MDP(8)) * g(3);];
taugX  = t1;
