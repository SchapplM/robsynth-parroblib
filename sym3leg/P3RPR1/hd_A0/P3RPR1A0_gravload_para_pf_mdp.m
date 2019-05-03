% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPR1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPR1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1A0_gravload_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1A0_gravload_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RPR1A0_gravload_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:18
% EndTime: 2019-05-03 14:58:19
% DurationCPUTime: 0.93s
% Computational Cost: add. (495->125), mult. (881->233), div. (54->3), fcn. (748->14), ass. (0->98)
t1328 = -MDP(3) + MDP(5);
t1323 = MDP(2) + MDP(4);
t1324 = pkin(1) * g(2);
t1296 = sin(qJ(1,3));
t1299 = cos(qJ(1,3));
t1308 = koppelP(3,2);
t1311 = koppelP(3,1);
t1261 = t1296 * t1311 - t1299 * t1308;
t1262 = t1296 * t1308 + t1299 * t1311;
t1293 = legFrame(3,3);
t1285 = sin(t1293);
t1288 = cos(t1293);
t1304 = xP(3);
t1291 = sin(t1304);
t1292 = cos(t1304);
t1223 = (t1261 * t1292 - t1262 * t1291) * t1288 + t1285 * (t1261 * t1291 + t1262 * t1292);
t1305 = 0.1e1 / qJ(2,3);
t1322 = t1223 * t1305;
t1297 = sin(qJ(1,2));
t1300 = cos(qJ(1,2));
t1309 = koppelP(2,2);
t1312 = koppelP(2,1);
t1263 = t1297 * t1312 - t1300 * t1309;
t1264 = t1297 * t1309 + t1300 * t1312;
t1294 = legFrame(2,3);
t1286 = sin(t1294);
t1289 = cos(t1294);
t1224 = (t1263 * t1292 - t1264 * t1291) * t1289 + t1286 * (t1263 * t1291 + t1264 * t1292);
t1306 = 0.1e1 / qJ(2,2);
t1321 = t1224 * t1306;
t1298 = sin(qJ(1,1));
t1301 = cos(qJ(1,1));
t1310 = koppelP(1,2);
t1313 = koppelP(1,1);
t1265 = t1298 * t1313 - t1301 * t1310;
t1266 = t1298 * t1310 + t1301 * t1313;
t1295 = legFrame(1,3);
t1287 = sin(t1295);
t1290 = cos(t1295);
t1225 = (t1265 * t1292 - t1266 * t1291) * t1290 + t1287 * (t1265 * t1291 + t1266 * t1292);
t1307 = 0.1e1 / qJ(2,1);
t1320 = t1225 * t1307;
t1247 = t1285 * t1299 + t1288 * t1296;
t1319 = t1247 * t1305;
t1248 = -t1285 * t1296 + t1288 * t1299;
t1318 = t1248 * t1305;
t1249 = t1286 * t1300 + t1289 * t1297;
t1317 = t1249 * t1306;
t1250 = -t1286 * t1297 + t1289 * t1300;
t1316 = t1250 * t1306;
t1251 = t1287 * t1301 + t1290 * t1298;
t1315 = t1251 * t1307;
t1252 = -t1287 * t1298 + t1290 * t1301;
t1314 = t1252 * t1307;
t1253 = g(1) * t1285 - g(2) * t1288;
t1256 = g(1) * t1288 + g(2) * t1285;
t1230 = t1253 * t1299 + t1256 * t1296;
t1229 = t1253 * t1296 - t1256 * t1299;
t1254 = g(1) * t1286 - g(2) * t1289;
t1257 = g(1) * t1289 + g(2) * t1286;
t1232 = t1254 * t1300 + t1257 * t1297;
t1231 = t1254 * t1297 - t1257 * t1300;
t1255 = g(1) * t1287 - g(2) * t1290;
t1258 = g(1) * t1290 + g(2) * t1287;
t1234 = t1255 * t1301 + t1258 * t1298;
t1233 = t1255 * t1298 - t1258 * t1301;
t1303 = pkin(1) * g(1);
t1302 = pkin(1) + pkin(2);
t1284 = g(1) * qJ(2,1) + t1324;
t1283 = -g(2) * qJ(2,1) + t1303;
t1282 = g(1) * qJ(2,2) + t1324;
t1281 = -g(2) * qJ(2,2) + t1303;
t1280 = g(1) * qJ(2,3) + t1324;
t1279 = -g(2) * qJ(2,3) + t1303;
t1278 = qJ(2,1) * t1313 + t1302 * t1310;
t1277 = qJ(2,2) * t1312 + t1302 * t1309;
t1276 = qJ(2,3) * t1311 + t1302 * t1308;
t1275 = -qJ(2,1) * t1310 + t1302 * t1313;
t1274 = -qJ(2,2) * t1309 + t1302 * t1312;
t1273 = -qJ(2,3) * t1308 + t1302 * t1311;
t1272 = qJ(2,1) * t1298 + t1301 * t1302;
t1271 = qJ(2,2) * t1297 + t1300 * t1302;
t1270 = qJ(2,3) * t1296 + t1299 * t1302;
t1269 = -qJ(2,1) * t1301 + t1298 * t1302;
t1268 = -qJ(2,2) * t1300 + t1297 * t1302;
t1267 = -qJ(2,3) * t1299 + t1296 * t1302;
t1260 = g(1) * t1292 + g(2) * t1291;
t1259 = g(1) * t1291 - g(2) * t1292;
t1246 = t1275 * t1298 - t1278 * t1301;
t1245 = t1274 * t1297 - t1277 * t1300;
t1244 = t1273 * t1296 - t1276 * t1299;
t1243 = t1275 * t1301 + t1278 * t1298;
t1242 = t1274 * t1300 + t1277 * t1297;
t1241 = t1273 * t1299 + t1276 * t1296;
t1228 = (t1283 * t1298 - t1284 * t1301) * t1290 + (t1283 * t1301 + t1284 * t1298) * t1287;
t1227 = (t1281 * t1297 - t1282 * t1300) * t1289 + (t1281 * t1300 + t1282 * t1297) * t1286;
t1226 = (t1279 * t1296 - t1280 * t1299) * t1288 + (t1279 * t1299 + t1280 * t1296) * t1285;
t1 = [((t1252 * t1228 - (-t1269 * t1287 + t1272 * t1290) * t1234) * t1307 + (t1250 * t1227 - (-t1268 * t1286 + t1271 * t1289) * t1232) * t1306 + (t1248 * t1226 - (-t1267 * t1285 + t1270 * t1288) * t1230) * t1305) * MDP(6) + (-t1259 * t1291 - t1260 * t1292) * MDP(10) + t1323 * (t1230 * t1318 + t1232 * t1316 + t1234 * t1314) + t1328 * (t1229 * t1318 + t1231 * t1316 + t1233 * t1314); ((t1251 * t1228 - (t1269 * t1290 + t1272 * t1287) * t1234) * t1307 + (t1249 * t1227 - (t1268 * t1289 + t1271 * t1286) * t1232) * t1306 + (t1247 * t1226 - (t1267 * t1288 + t1270 * t1285) * t1230) * t1305) * MDP(6) + (t1259 * t1292 - t1260 * t1291) * MDP(10) + t1323 * (t1230 * t1319 + t1232 * t1317 + t1234 * t1315) + t1328 * (t1229 * t1319 + t1231 * t1317 + t1233 * t1315); ((t1225 * t1228 - ((-t1243 * t1291 + t1246 * t1292) * t1290 + (t1243 * t1292 + t1246 * t1291) * t1287) * t1234) * t1307 + (t1224 * t1227 - ((-t1242 * t1291 + t1245 * t1292) * t1289 + (t1242 * t1292 + t1245 * t1291) * t1286) * t1232) * t1306 + (t1223 * t1226 - ((-t1241 * t1291 + t1244 * t1292) * t1288 + (t1241 * t1292 + t1244 * t1291) * t1285) * t1230) * t1305) * MDP(6) + t1259 * MDP(8) + t1260 * MDP(9) + t1323 * (t1230 * t1322 + t1232 * t1321 + t1234 * t1320) + t1328 * (t1229 * t1322 + t1231 * t1321 + t1233 * t1320);];
taugX  = t1;
