% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR2G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G2A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:45
% EndTime: 2020-03-09 21:21:47
% DurationCPUTime: 1.89s
% Computational Cost: add. (9465->180), mult. (10029->365), div. (2034->9), fcn. (7461->33), ass. (0->203)
t1316 = sin(qJ(3,3));
t1305 = 0.1e1 / t1316 ^ 2;
t1322 = cos(qJ(3,3));
t1445 = t1305 * t1322;
t1318 = sin(qJ(3,2));
t1308 = 0.1e1 / t1318 ^ 2;
t1324 = cos(qJ(3,2));
t1444 = t1308 * t1324;
t1320 = sin(qJ(3,1));
t1311 = 0.1e1 / t1320 ^ 2;
t1326 = cos(qJ(3,1));
t1443 = t1311 * t1326;
t1332 = 0.1e1 / pkin(2);
t1442 = (t1322 * pkin(1) + pkin(2)) * t1332;
t1441 = (t1324 * pkin(1) + pkin(2)) * t1332;
t1440 = (t1326 * pkin(1) + pkin(2)) * t1332;
t1304 = 0.1e1 / t1316;
t1439 = t1304 * t1445;
t1307 = 0.1e1 / t1318;
t1438 = t1307 * t1444;
t1310 = 0.1e1 / t1320;
t1437 = t1310 * t1443;
t1328 = xDP(3);
t1436 = -2 * t1328;
t1313 = legFrame(3,2);
t1426 = qJ(2,3) + qJ(3,3);
t1283 = t1313 + t1426;
t1284 = -t1313 + t1426;
t1295 = sin(t1426);
t1329 = xDP(2);
t1330 = xDP(1);
t1268 = t1295 * t1436 + (sin(t1283) - sin(t1284)) * t1330 + (cos(t1284) + cos(t1283)) * t1329;
t1435 = -t1268 / 0.2e1;
t1314 = legFrame(2,2);
t1427 = qJ(2,2) + qJ(3,2);
t1285 = t1314 + t1427;
t1286 = -t1314 + t1427;
t1296 = sin(t1427);
t1269 = t1296 * t1436 + (-sin(t1286) + sin(t1285)) * t1330 + (cos(t1286) + cos(t1285)) * t1329;
t1434 = -t1269 / 0.2e1;
t1315 = legFrame(1,2);
t1428 = qJ(2,1) + qJ(3,1);
t1287 = t1315 + t1428;
t1288 = -t1315 + t1428;
t1297 = sin(t1428);
t1270 = t1297 * t1436 + (sin(t1287) - sin(t1288)) * t1330 + (cos(t1288) + cos(t1287)) * t1329;
t1433 = -t1270 / 0.2e1;
t1432 = pkin(2) * t1322;
t1431 = pkin(2) * t1324;
t1430 = pkin(2) * t1326;
t1429 = pkin(2) * t1328;
t1334 = 0.1e1 / pkin(1) ^ 2;
t1400 = t1305 * t1334;
t1333 = 0.1e1 / pkin(1);
t1392 = t1332 * t1333;
t1298 = sin(t1313);
t1301 = cos(t1313);
t1274 = t1298 * t1330 + t1301 * t1329;
t1292 = pkin(1) + t1432;
t1317 = sin(qJ(2,3));
t1323 = cos(qJ(2,3));
t1259 = (-t1274 * t1292 + t1316 * t1429) * t1323 + t1317 * (pkin(2) * t1274 * t1316 + t1292 * t1328);
t1419 = t1259 * t1304;
t1256 = t1392 * t1419;
t1416 = t1268 * t1304;
t1367 = t1416 / 0.2e1;
t1262 = t1333 * t1367;
t1250 = t1262 + t1256;
t1422 = t1250 * t1259;
t1244 = t1400 * t1422;
t1401 = t1304 * t1322;
t1247 = -pkin(2) * t1250 + t1401 * t1435;
t1352 = t1400 * t1435;
t1235 = t1247 * t1352 + t1244;
t1331 = pkin(2) ^ 2;
t1238 = (t1250 * t1331 + (0.2e1 * (t1262 + t1256 / 0.2e1) * t1432 + t1367) * pkin(1)) * t1332 * t1352;
t1229 = -t1244 * t1442 + t1235 + t1238;
t1425 = t1229 * t1304;
t1398 = t1308 * t1334;
t1299 = sin(t1314);
t1302 = cos(t1314);
t1275 = t1299 * t1330 + t1302 * t1329;
t1293 = pkin(1) + t1431;
t1319 = sin(qJ(2,2));
t1325 = cos(qJ(2,2));
t1260 = (-t1275 * t1293 + t1318 * t1429) * t1325 + t1319 * (pkin(2) * t1275 * t1318 + t1293 * t1328);
t1418 = t1260 * t1307;
t1257 = t1392 * t1418;
t1415 = t1269 * t1307;
t1366 = t1415 / 0.2e1;
t1263 = t1333 * t1366;
t1251 = t1263 + t1257;
t1421 = t1251 * t1260;
t1245 = t1398 * t1421;
t1399 = t1307 * t1324;
t1248 = -pkin(2) * t1251 + t1399 * t1434;
t1351 = t1398 * t1434;
t1236 = t1248 * t1351 + t1245;
t1239 = (t1251 * t1331 + (0.2e1 * (t1263 + t1257 / 0.2e1) * t1431 + t1366) * pkin(1)) * t1332 * t1351;
t1230 = -t1245 * t1441 + t1236 + t1239;
t1424 = t1230 * t1307;
t1396 = t1311 * t1334;
t1300 = sin(t1315);
t1303 = cos(t1315);
t1276 = t1300 * t1330 + t1303 * t1329;
t1294 = pkin(1) + t1430;
t1321 = sin(qJ(2,1));
t1327 = cos(qJ(2,1));
t1261 = (-t1276 * t1294 + t1320 * t1429) * t1327 + t1321 * (pkin(2) * t1276 * t1320 + t1294 * t1328);
t1417 = t1261 * t1310;
t1258 = t1392 * t1417;
t1414 = t1270 * t1310;
t1365 = t1414 / 0.2e1;
t1264 = t1333 * t1365;
t1252 = t1264 + t1258;
t1420 = t1252 * t1261;
t1246 = t1396 * t1420;
t1397 = t1310 * t1326;
t1249 = -t1252 * pkin(2) + t1397 * t1433;
t1350 = t1396 * t1433;
t1237 = t1249 * t1350 + t1246;
t1240 = (t1252 * t1331 + (0.2e1 * (t1264 + t1258 / 0.2e1) * t1430 + t1365) * pkin(1)) * t1332 * t1350;
t1231 = -t1246 * t1440 + t1237 + t1240;
t1423 = t1231 * t1310;
t1395 = t1317 * t1316;
t1271 = -pkin(2) * t1395 + t1292 * t1323;
t1413 = t1271 * t1235;
t1394 = t1319 * t1318;
t1272 = -pkin(2) * t1394 + t1293 * t1325;
t1412 = t1272 * t1236;
t1393 = t1321 * t1320;
t1273 = -pkin(2) * t1393 + t1294 * t1327;
t1411 = t1273 * t1237;
t1410 = t1295 * t1304;
t1409 = t1296 * t1307;
t1408 = t1297 * t1310;
t1277 = t1323 * t1322 - t1395;
t1407 = t1298 * t1277;
t1278 = t1325 * t1324 - t1394;
t1406 = t1299 * t1278;
t1279 = t1327 * t1326 - t1393;
t1405 = t1300 * t1279;
t1404 = t1301 * t1277;
t1403 = t1302 * t1278;
t1402 = t1303 * t1279;
t1391 = t1271 * t1425;
t1390 = t1272 * t1424;
t1389 = t1273 * t1423;
t1253 = t1333 * t1416 + t1256;
t1388 = t1253 * t1419;
t1254 = t1333 * t1415 + t1257;
t1387 = t1254 * t1418;
t1255 = t1333 * t1414 + t1258;
t1386 = t1255 * t1417;
t1385 = t1304 * t1407;
t1384 = t1304 * t1404;
t1383 = t1307 * t1406;
t1382 = t1307 * t1403;
t1381 = t1310 * t1405;
t1380 = t1310 * t1402;
t1232 = t1238 + 0.2e1 * t1244 + (-t1247 * t1268 - t1422 * t1442) * t1400;
t1379 = t1232 * t1401;
t1378 = t1235 * t1401;
t1233 = t1239 + 0.2e1 * t1245 + (-t1248 * t1269 - t1421 * t1441) * t1398;
t1377 = t1233 * t1399;
t1376 = t1236 * t1399;
t1234 = t1240 + 0.2e1 * t1246 + (-t1249 * t1270 - t1420 * t1440) * t1396;
t1375 = t1234 * t1397;
t1374 = t1237 * t1397;
t1265 = t1268 ^ 2;
t1373 = -t1265 * t1271 / 0.4e1;
t1280 = t1317 * pkin(1) + pkin(2) * t1295;
t1372 = t1265 * t1280 / 0.4e1;
t1266 = t1269 ^ 2;
t1371 = -t1266 * t1272 / 0.4e1;
t1281 = t1319 * pkin(1) + pkin(2) * t1296;
t1370 = t1266 * t1281 / 0.4e1;
t1267 = t1270 ^ 2;
t1369 = -t1267 * t1273 / 0.4e1;
t1282 = t1321 * pkin(1) + pkin(2) * t1297;
t1368 = t1267 * t1282 / 0.4e1;
t1364 = t1277 * t1388;
t1363 = t1253 * t1259 * t1445;
t1362 = t1278 * t1387;
t1361 = t1254 * t1260 * t1444;
t1360 = t1279 * t1386;
t1359 = t1255 * t1261 * t1443;
t1358 = t1271 * t1378;
t1357 = t1272 * t1376;
t1356 = t1273 * t1374;
t1355 = t1277 * t1379;
t1354 = t1278 * t1377;
t1353 = t1279 * t1375;
t1349 = t1305 * t1373;
t1348 = t1308 * t1371;
t1347 = t1311 * t1369;
t1346 = t1277 * t1363;
t1345 = t1278 * t1361;
t1344 = t1279 * t1359;
t1343 = t1373 * t1439;
t1342 = t1371 * t1438;
t1341 = t1369 * t1437;
t1 = [(t1298 * t1355 + t1299 * t1354 + t1300 * t1353) * MDP(6) + (-t1232 * t1407 - t1233 * t1406 - t1234 * t1405) * MDP(7) + ((t1235 * t1385 + t1236 * t1383 + t1237 * t1381) * MDP(2) + (t1229 * t1385 + t1230 * t1383 + t1231 * t1381) * MDP(5)) * t1333 + ((-t1298 * t1358 - t1299 * t1357 - t1300 * t1356) * MDP(6) + (t1298 * t1413 + t1299 * t1412 + t1300 * t1411) * MDP(7) + ((t1298 * t1349 + t1299 * t1348 + t1300 * t1347) * MDP(6) + (t1298 * t1343 + t1299 * t1342 + t1300 * t1341) * MDP(7)) * t1334 + ((-t1298 * t1391 - t1299 * t1390 - t1300 * t1389) * MDP(5) + (-t1298 * t1364 - t1299 * t1362 - t1300 * t1360) * MDP(6) + (-t1298 * t1346 - t1299 * t1345 - t1300 * t1344) * MDP(7)) * t1333) * t1332; (t1301 * t1355 + t1302 * t1354 + t1303 * t1353) * MDP(6) + (-t1232 * t1404 - t1233 * t1403 - t1234 * t1402) * MDP(7) + ((t1235 * t1384 + t1236 * t1382 + t1237 * t1380) * MDP(2) + (t1229 * t1384 + t1230 * t1382 + t1231 * t1380) * MDP(5)) * t1333 + ((-t1301 * t1358 - t1302 * t1357 - t1303 * t1356) * MDP(6) + (t1301 * t1413 + t1302 * t1412 + t1303 * t1411) * MDP(7) + ((t1301 * t1349 + t1302 * t1348 + t1303 * t1347) * MDP(6) + (t1301 * t1343 + t1302 * t1342 + t1303 * t1341) * MDP(7)) * t1334 + ((-t1301 * t1391 - t1302 * t1390 - t1303 * t1389) * MDP(5) + (-t1301 * t1364 - t1302 * t1362 - t1303 * t1360) * MDP(6) + (-t1301 * t1346 - t1302 * t1345 - t1303 * t1344) * MDP(7)) * t1333) * t1332; (-t1295 * t1379 - t1296 * t1377 - t1297 * t1375) * MDP(6) + (t1295 * t1232 + t1296 * t1233 + t1297 * t1234) * MDP(7) + ((-t1235 * t1410 - t1236 * t1409 - t1237 * t1408) * MDP(2) + (-t1229 * t1410 - t1230 * t1409 - t1231 * t1408) * MDP(5)) * t1333 + ((t1280 * t1378 + t1281 * t1376 + t1282 * t1374) * MDP(6) + (-t1280 * t1235 - t1281 * t1236 - t1282 * t1237) * MDP(7) + ((t1305 * t1372 + t1308 * t1370 + t1311 * t1368) * MDP(6) + (t1368 * t1437 + t1370 * t1438 + t1372 * t1439) * MDP(7)) * t1334 + ((t1280 * t1425 + t1281 * t1424 + t1282 * t1423) * MDP(5) + (t1295 * t1388 + t1296 * t1387 + t1297 * t1386) * MDP(6) + (t1295 * t1363 + t1296 * t1361 + t1297 * t1359) * MDP(7)) * t1333) * t1332;];
taucX  = t1;
