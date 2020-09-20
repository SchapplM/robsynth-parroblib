% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR1G3A0
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
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G3A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:33
% EndTime: 2020-03-09 21:02:35
% DurationCPUTime: 2.53s
% Computational Cost: add. (1254->210), mult. (4817->471), div. (1707->17), fcn. (4095->18), ass. (0->207)
t1305 = sin(qJ(2,3));
t1269 = 0.1e1 / t1305;
t1310 = cos(qJ(3,3));
t1333 = t1310 ^ 2;
t1283 = 0.1e1 / t1333;
t1304 = sin(qJ(3,3));
t1268 = t1304 ^ 2;
t1285 = 0.1e1 / t1333 ^ 2;
t1449 = t1268 * t1285;
t1469 = t1269 * (t1283 + t1449);
t1307 = sin(qJ(2,2));
t1274 = 0.1e1 / t1307;
t1312 = cos(qJ(3,2));
t1337 = t1312 ^ 2;
t1289 = 0.1e1 / t1337;
t1306 = sin(qJ(3,2));
t1273 = t1306 ^ 2;
t1291 = 0.1e1 / t1337 ^ 2;
t1445 = t1273 * t1291;
t1468 = t1274 * (t1289 + t1445);
t1309 = sin(qJ(2,1));
t1279 = 0.1e1 / t1309;
t1314 = cos(qJ(3,1));
t1341 = t1314 ^ 2;
t1295 = 0.1e1 / t1341;
t1308 = sin(qJ(3,1));
t1278 = t1308 ^ 2;
t1297 = 0.1e1 / t1341 ^ 2;
t1441 = t1278 * t1297;
t1467 = t1279 * (t1295 + t1441);
t1319 = 0.1e1 / pkin(2);
t1466 = 2 * MDP(6);
t1282 = 0.1e1 / t1310;
t1288 = 0.1e1 / t1312;
t1294 = 0.1e1 / t1314;
t1320 = 0.1e1 / pkin(2) ^ 2;
t1270 = 0.1e1 / t1305 ^ 2;
t1275 = 0.1e1 / t1307 ^ 2;
t1280 = 0.1e1 / t1309 ^ 2;
t1316 = xDP(3);
t1300 = t1316 ^ 2;
t1465 = MDP(7) * t1300;
t1301 = legFrame(3,2);
t1261 = sin(t1301);
t1264 = cos(t1301);
t1317 = xDP(2);
t1318 = xDP(1);
t1311 = cos(qJ(2,3));
t1420 = t1311 * t1316;
t1249 = -t1304 * t1420 + (t1261 * t1317 - t1264 * t1318) * t1310;
t1284 = t1282 * t1283;
t1434 = t1285 * t1320;
t1382 = t1270 * t1434;
t1446 = t1269 * t1311;
t1389 = t1249 * t1446;
t1417 = t1316 * t1319 ^ 2;
t1423 = t1304 * t1305;
t1448 = t1269 * t1282;
t1228 = -(-t1316 * t1423 + t1389) * t1249 * t1382 - (-t1249 * t1304 + t1420) * t1284 * t1417 * t1448;
t1464 = t1228 * t1311;
t1302 = legFrame(2,2);
t1262 = sin(t1302);
t1265 = cos(t1302);
t1313 = cos(qJ(2,2));
t1419 = t1313 * t1316;
t1250 = -t1306 * t1419 + (t1262 * t1317 - t1265 * t1318) * t1312;
t1290 = t1288 * t1289;
t1430 = t1291 * t1320;
t1380 = t1275 * t1430;
t1442 = t1274 * t1313;
t1387 = t1250 * t1442;
t1422 = t1306 * t1307;
t1444 = t1274 * t1288;
t1229 = -(-t1316 * t1422 + t1387) * t1250 * t1380 - (-t1250 * t1306 + t1419) * t1290 * t1417 * t1444;
t1463 = t1229 * t1313;
t1303 = legFrame(1,2);
t1263 = sin(t1303);
t1266 = cos(t1303);
t1315 = cos(qJ(2,1));
t1418 = t1315 * t1316;
t1251 = -t1308 * t1418 + (t1263 * t1317 - t1266 * t1318) * t1314;
t1296 = t1294 * t1295;
t1426 = t1297 * t1320;
t1378 = t1280 * t1426;
t1438 = t1279 * t1315;
t1385 = t1251 * t1438;
t1421 = t1308 * t1309;
t1440 = t1279 * t1294;
t1230 = -(-t1316 * t1421 + t1385) * t1251 * t1378 - (-t1251 * t1308 + t1418) * t1296 * t1417 * t1440;
t1462 = t1230 * t1315;
t1246 = t1249 ^ 2;
t1271 = t1269 * t1270;
t1458 = t1246 * t1271;
t1237 = (t1269 * t1300 + t1458) * t1319 * t1284;
t1461 = t1237 * t1282;
t1247 = t1250 ^ 2;
t1276 = t1274 * t1275;
t1456 = t1247 * t1276;
t1238 = (t1274 * t1300 + t1456) * t1319 * t1290;
t1460 = t1238 * t1288;
t1248 = t1251 ^ 2;
t1281 = t1279 * t1280;
t1454 = t1248 * t1281;
t1239 = (t1279 * t1300 + t1454) * t1319 * t1296;
t1459 = t1239 * t1294;
t1457 = t1246 * t1320;
t1455 = t1247 * t1320;
t1453 = t1248 * t1320;
t1452 = t1249 * t1270;
t1451 = t1250 * t1275;
t1450 = t1251 * t1280;
t1287 = t1311 ^ 2;
t1447 = t1269 * t1287;
t1293 = t1313 ^ 2;
t1443 = t1274 * t1293;
t1299 = t1315 ^ 2;
t1439 = t1279 * t1299;
t1437 = t1282 * t1304;
t1436 = t1284 * t1304;
t1435 = t1285 * t1304;
t1433 = t1288 * t1306;
t1432 = t1290 * t1306;
t1431 = t1291 * t1306;
t1429 = t1294 * t1308;
t1428 = t1296 * t1308;
t1427 = t1297 * t1308;
t1425 = t1300 * t1320;
t1243 = t1246 * t1382;
t1240 = t1283 * t1425 + t1243;
t1408 = -0.2e1 * t1316 * t1320;
t1350 = t1389 * t1408;
t1377 = t1305 * t1425;
t1416 = -t1268 * t1284 * t1377 + t1350 * t1436 + (-t1240 * t1305 + t1464) * t1310;
t1244 = t1247 * t1380;
t1241 = t1289 * t1425 + t1244;
t1349 = t1387 * t1408;
t1376 = t1307 * t1425;
t1415 = -t1273 * t1290 * t1376 + t1349 * t1432 + (-t1241 * t1307 + t1463) * t1312;
t1245 = t1248 * t1378;
t1242 = t1295 * t1425 + t1245;
t1348 = t1385 * t1408;
t1375 = t1309 * t1425;
t1414 = -t1278 * t1296 * t1375 + t1348 * t1428 + (-t1242 * t1309 + t1462) * t1314;
t1413 = (-t1283 * t1377 - t1464) * t1304 + t1240 * t1423 + t1283 * t1350;
t1412 = (-t1289 * t1376 - t1463) * t1306 + t1241 * t1422 + t1289 * t1349;
t1411 = (-t1295 * t1375 - t1462) * t1308 + t1242 * t1421 + t1295 * t1348;
t1410 = 2 * MDP(5);
t1409 = 0.2e1 * t1316;
t1407 = t1228 * t1448;
t1406 = t1228 * t1269 * t1304;
t1405 = t1228 * t1446;
t1404 = t1229 * t1444;
t1403 = t1229 * t1274 * t1306;
t1402 = t1229 * t1442;
t1401 = t1230 * t1440;
t1400 = t1230 * t1279 * t1308;
t1399 = t1230 * t1438;
t1398 = t1237 * t1446;
t1397 = t1237 * t1283 * t1319;
t1396 = t1238 * t1442;
t1395 = t1238 * t1289 * t1319;
t1394 = t1239 * t1438;
t1393 = t1239 * t1295 * t1319;
t1286 = t1282 * t1285;
t1392 = t1286 * t1457;
t1292 = t1288 * t1291;
t1391 = t1292 * t1455;
t1298 = t1294 * t1297;
t1390 = t1298 * t1453;
t1388 = t1311 * t1452;
t1386 = t1313 * t1451;
t1384 = t1315 * t1450;
t1383 = t1269 * t1437;
t1381 = t1274 * t1433;
t1379 = t1279 * t1429;
t1371 = t1268 * t1407;
t1370 = t1283 * t1405;
t1369 = t1273 * t1404;
t1368 = t1289 * t1402;
t1367 = t1278 * t1401;
t1366 = t1295 * t1399;
t1365 = t1282 * t1398;
t1364 = t1288 * t1396;
t1363 = t1294 * t1394;
t1258 = -0.1e1 + 0.2e1 * t1333;
t1362 = t1258 * t1285 * t1452;
t1361 = t1436 * t1452;
t1259 = -0.1e1 + 0.2e1 * t1337;
t1360 = t1259 * t1291 * t1451;
t1359 = t1432 * t1451;
t1260 = -0.1e1 + 0.2e1 * t1341;
t1358 = t1260 * t1297 * t1450;
t1357 = t1428 * t1450;
t1356 = t1261 * t1365;
t1355 = t1262 * t1364;
t1354 = t1263 * t1363;
t1353 = t1264 * t1365;
t1352 = t1265 * t1364;
t1351 = t1266 * t1363;
t1347 = -MDP(3) * t1243 + (-t1311 * t1434 * t1458 - t1228) * MDP(4) + (t1237 * MDP(1) + t1416 * MDP(10) + t1413 * MDP(11) + MDP(3) * t1464) * t1269;
t1346 = -MDP(3) * t1244 + (-t1313 * t1430 * t1456 - t1229) * MDP(4) + (t1238 * MDP(1) + t1415 * MDP(10) + t1412 * MDP(11) + MDP(3) * t1463) * t1274;
t1345 = -MDP(3) * t1245 + (-t1315 * t1426 * t1454 - t1230) * MDP(4) + (t1239 * MDP(1) + t1414 * MDP(10) + t1411 * MDP(11) + MDP(3) * t1462) * t1279;
t1321 = t1319 * t1320;
t1277 = t1308 * t1278;
t1272 = t1306 * t1273;
t1267 = t1304 * t1268;
t1 = [t1345 * (t1263 * t1309 + t1266 * t1315) + t1346 * (t1262 * t1307 + t1265 * t1313) + t1347 * (t1261 * t1305 + t1264 * t1311) + (((-t1264 * t1361 - t1265 * t1359 - t1266 * t1357) * t1410 + (-t1264 * t1362 - t1265 * t1360 - t1266 * t1358) * t1466) * t1316 + (-t1264 * t1469 - t1265 * t1468 - t1266 * t1467) * t1465) * t1321 + ((-t1264 * t1407 - t1265 * t1404 - t1266 * t1401) * MDP(2) + (-t1351 - t1352 - t1353) * MDP(3) + (t1264 * t1461 + t1265 * t1460 + t1266 * t1459) * MDP(4) + (-t1264 * t1371 - t1265 * t1369 - t1266 * t1367) * MDP(5) + (-t1264 * t1406 - t1265 * t1403 - t1266 * t1400) * t1466 + (-t1264 * t1398 - t1265 * t1396 - t1266 * t1394) * MDP(10) + (t1304 * t1353 + t1306 * t1352 + t1308 * t1351) * MDP(11)) * t1319; t1345 * (-t1263 * t1315 + t1266 * t1309) + t1346 * (-t1262 * t1313 + t1265 * t1307) + t1347 * (-t1261 * t1311 + t1264 * t1305) + (((t1261 * t1361 + t1262 * t1359 + t1263 * t1357) * t1410 + (t1261 * t1362 + t1262 * t1360 + t1263 * t1358) * t1466) * t1316 + (t1261 * t1469 + t1262 * t1468 + t1263 * t1467) * t1465) * t1321 + ((t1261 * t1407 + t1262 * t1404 + t1263 * t1401) * MDP(2) + (t1354 + t1355 + t1356) * MDP(3) + (-t1261 * t1461 - t1262 * t1460 - t1263 * t1459) * MDP(4) + (t1261 * t1371 + t1262 * t1369 + t1263 * t1367) * MDP(5) + (t1261 * t1406 + t1262 * t1403 + t1263 * t1400) * t1466 + (t1261 * t1398 + t1262 * t1396 + t1263 * t1394) * MDP(10) + (-t1304 * t1356 - t1306 * t1355 - t1308 * t1354) * MDP(11)) * t1319; (t1237 * t1383 + t1238 * t1381 + t1239 * t1379) * MDP(1) + ((-t1280 * t1390 + (t1294 * t1462 - t1299 * t1393) * t1279) * t1308 + (-t1275 * t1391 + (t1288 * t1463 - t1293 * t1395) * t1274) * t1306 + (-t1270 * t1392 + (t1282 * t1464 - t1287 * t1397) * t1269) * t1304) * MDP(3) + ((-t1230 * t1294 + (-t1281 * t1390 + t1393) * t1315) * t1308 + (-t1229 * t1288 + (-t1276 * t1391 + t1395) * t1313) * t1306 + (-t1228 * t1282 + (-t1271 * t1392 + t1397) * t1311) * t1304) * MDP(4) + (t1414 * t1379 + t1415 * t1381 + t1416 * t1383) * MDP(10) + (t1411 * t1379 + t1412 * t1381 + t1413 * t1383) * MDP(11) + ((-t1304 * t1370 - t1306 * t1368 - t1308 * t1366) * MDP(2) + (-t1267 * t1370 - t1272 * t1368 - t1277 * t1366) * MDP(5) + ((-0.2e1 * t1280 * t1295 * t1453 - 0.2e1 * t1278 * t1399 + t1245) * t1294 + (-0.2e1 * t1275 * t1289 * t1455 - 0.2e1 * t1273 * t1402 + t1244) * t1288 + (-0.2e1 * t1270 * t1283 * t1457 - 0.2e1 * t1268 * t1405 + t1243) * t1282) * MDP(6) + (t1228 * t1437 + t1229 * t1433 + t1230 * t1429) * MDP(7) + (t1228 + t1229 + t1230) * MDP(8) + ((-t1309 - t1439) * t1239 * t1429 + (-t1307 - t1443) * t1238 * t1433 + (-t1305 - t1447) * t1237 * t1437) * MDP(10) + ((t1278 * t1295 * t1439 - t1309) * t1239 + (t1273 * t1289 * t1443 - t1307) * t1238 + (t1268 * t1283 * t1447 - t1305) * t1237) * MDP(11)) * t1319 + ((((-t1277 * t1298 - t1428) * t1438 + (-t1272 * t1292 - t1432) * t1442 + (-t1267 * t1286 - t1436) * t1446) * MDP(7) + (t1427 + t1431 + t1435) * MDP(9)) * t1300 + (-t1246 * t1270 * t1435 - t1247 * t1275 * t1431 - t1248 * t1280 * t1427 + (-t1384 * t1441 - t1386 * t1445 - t1388 * t1449) * t1409) * MDP(5) + (-t1258 * t1286 * t1304 * t1388 - t1259 * t1292 * t1306 * t1386 - t1260 * t1298 * t1308 * t1384) * t1409 * MDP(6)) * t1321;];
taucX  = t1;
