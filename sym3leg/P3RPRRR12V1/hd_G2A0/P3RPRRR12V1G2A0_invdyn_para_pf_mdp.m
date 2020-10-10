% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR12V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR12V1G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:55
% EndTime: 2020-08-06 18:24:59
% DurationCPUTime: 4.14s
% Computational Cost: add. (15597->344), mult. (25502->606), div. (3228->14), fcn. (18972->18), ass. (0->255)
t1389 = sin(qJ(3,3));
t1374 = 0.1e1 / t1389 ^ 2;
t1383 = legFrame(3,2);
t1356 = sin(t1383);
t1359 = cos(t1383);
t1402 = xDP(2);
t1403 = xDP(1);
t1321 = t1356 * t1403 + t1359 * t1402;
t1413 = 0.1e1 / pkin(3) ^ 2;
t1520 = t1321 ^ 2 * t1413;
t1472 = t1374 * t1520;
t1391 = sin(qJ(3,2));
t1376 = 0.1e1 / t1391 ^ 2;
t1384 = legFrame(2,2);
t1357 = sin(t1384);
t1360 = cos(t1384);
t1323 = t1357 * t1403 + t1360 * t1402;
t1519 = t1323 ^ 2 * t1413;
t1470 = t1376 * t1519;
t1393 = sin(qJ(3,1));
t1378 = 0.1e1 / t1393 ^ 2;
t1385 = legFrame(1,2);
t1358 = sin(t1385);
t1361 = cos(t1385);
t1325 = t1358 * t1403 + t1361 * t1402;
t1518 = t1325 ^ 2 * t1413;
t1468 = t1378 * t1518;
t1551 = (-pkin(5) - pkin(6));
t1372 = pkin(1) - t1551;
t1401 = xDP(3);
t1348 = t1401 * t1372;
t1365 = pkin(3) * t1401;
t1390 = sin(qJ(1,3));
t1396 = cos(qJ(1,3));
t1454 = t1356 * t1402 - t1359 * t1403;
t1395 = cos(qJ(3,3));
t1544 = pkin(3) * t1396;
t1466 = (t1395 + 0.1e1) * (t1395 - 0.1e1) * t1544;
t1487 = t1395 * t1321;
t1379 = t1395 ^ 2;
t1564 = -t1379 + 0.1e1;
t1258 = ((t1454 * qJ(2,3) + t1348) * t1396 + (qJ(2,3) * t1401 - t1454 * t1372) * t1390 + pkin(3) * t1487) * t1389 - t1454 * t1466 + qJ(2,3) * t1487 + t1564 * t1390 * t1365;
t1303 = -t1454 * t1390 + t1396 * t1401;
t1333 = g(1) * t1359 - g(2) * t1356;
t1312 = g(3) * t1396 + t1333 * t1390;
t1547 = pkin(3) * t1389;
t1349 = qJ(2,3) + t1547;
t1340 = 0.1e1 / t1349 ^ 2;
t1373 = 0.1e1 / t1389;
t1510 = t1340 * t1373;
t1573 = 0.2e1 * t1303 * t1258 * t1510 - t1312;
t1392 = sin(qJ(1,2));
t1398 = cos(qJ(1,2));
t1453 = t1357 * t1402 - t1360 * t1403;
t1397 = cos(qJ(3,2));
t1543 = pkin(3) * t1398;
t1465 = (t1397 + 0.1e1) * (t1397 - 0.1e1) * t1543;
t1486 = t1397 * t1323;
t1380 = t1397 ^ 2;
t1563 = -t1380 + 0.1e1;
t1259 = ((t1453 * qJ(2,2) + t1348) * t1398 + (qJ(2,2) * t1401 - t1453 * t1372) * t1392 + pkin(3) * t1486) * t1391 - t1453 * t1465 + qJ(2,2) * t1486 + t1563 * t1392 * t1365;
t1304 = -t1453 * t1392 + t1398 * t1401;
t1334 = g(1) * t1360 - g(2) * t1357;
t1313 = g(3) * t1398 + t1334 * t1392;
t1546 = pkin(3) * t1391;
t1350 = qJ(2,2) + t1546;
t1343 = 0.1e1 / t1350 ^ 2;
t1375 = 0.1e1 / t1391;
t1508 = t1343 * t1375;
t1572 = 0.2e1 * t1304 * t1259 * t1508 - t1313;
t1394 = sin(qJ(1,1));
t1400 = cos(qJ(1,1));
t1452 = t1358 * t1402 - t1361 * t1403;
t1399 = cos(qJ(3,1));
t1542 = pkin(3) * t1400;
t1464 = (t1399 + 0.1e1) * (t1399 - 0.1e1) * t1542;
t1485 = t1399 * t1325;
t1381 = t1399 ^ 2;
t1562 = -t1381 + 0.1e1;
t1260 = ((t1452 * qJ(2,1) + t1348) * t1400 + (qJ(2,1) * t1401 - t1452 * t1372) * t1394 + pkin(3) * t1485) * t1393 - t1452 * t1464 + qJ(2,1) * t1485 + t1562 * t1394 * t1365;
t1305 = -t1452 * t1394 + t1400 * t1401;
t1335 = g(1) * t1361 - g(2) * t1358;
t1314 = g(3) * t1400 + t1335 * t1394;
t1545 = pkin(3) * t1393;
t1351 = qJ(2,1) + t1545;
t1346 = 0.1e1 / t1351 ^ 2;
t1377 = 0.1e1 / t1393;
t1506 = t1346 * t1377;
t1571 = 0.2e1 * t1305 * t1260 * t1506 - t1314;
t1327 = t1333 * t1396;
t1339 = 0.1e1 / t1349;
t1341 = t1339 * t1340;
t1386 = xDDP(3);
t1387 = xDDP(2);
t1388 = xDDP(1);
t1525 = t1303 * t1372;
t1526 = t1303 * t1340;
t1243 = (-t1258 * t1341 + 0.2e1 * t1340 * t1487) * t1373 * t1303 + (t1396 * t1386 - (t1258 * t1373 - t1525) * t1526 + (-t1356 * t1387 + t1359 * t1388) * t1390) * t1339;
t1550 = pkin(1) * t1243;
t1570 = t1327 - t1550;
t1328 = t1334 * t1398;
t1342 = 0.1e1 / t1350;
t1344 = t1342 * t1343;
t1523 = t1304 * t1372;
t1524 = t1304 * t1343;
t1244 = (-t1259 * t1344 + 0.2e1 * t1343 * t1486) * t1375 * t1304 + (t1398 * t1386 - (t1259 * t1375 - t1523) * t1524 + (-t1357 * t1387 + t1360 * t1388) * t1392) * t1342;
t1549 = pkin(1) * t1244;
t1569 = t1328 - t1549;
t1329 = t1335 * t1400;
t1345 = 0.1e1 / t1351;
t1347 = t1345 * t1346;
t1521 = t1305 * t1372;
t1522 = t1305 * t1346;
t1245 = (-t1260 * t1347 + 0.2e1 * t1346 * t1485) * t1377 * t1305 + (t1400 * t1386 - (t1260 * t1377 - t1521) * t1522 + (-t1358 * t1387 + t1361 * t1388) * t1394) * t1345;
t1548 = pkin(1) * t1245;
t1568 = t1329 - t1548;
t1404 = pkin(1) + pkin(5);
t1362 = g(3) * t1390;
t1337 = qJ(2,3) * t1396 - t1372 * t1390;
t1504 = t1356 * t1395;
t1289 = (pkin(3) * t1504 - t1337 * t1359) * t1389 + t1359 * t1466 + qJ(2,3) * t1504;
t1498 = t1359 * t1395;
t1292 = (pkin(3) * t1498 + t1337 * t1356) * t1389 + t1564 * t1356 * t1544 + qJ(2,3) * t1498;
t1315 = t1349 * t1390 + t1372 * t1396;
t1406 = qJ(2,3) ^ 2;
t1411 = pkin(3) ^ 2;
t1445 = -(pkin(6) ^ 2) - t1411 + (-(2 * pkin(6)) - pkin(5)) * pkin(5) + ((2 * t1551) - pkin(1)) * pkin(1);
t1493 = t1372 * t1373;
t1481 = t1258 * t1493;
t1412 = 0.1e1 / pkin(3);
t1490 = t1373 * t1412;
t1517 = t1321 * t1339;
t1448 = (-t1487 * t1493 + (t1481 + (-0.2e1 * qJ(2,3) * t1547 + t1379 * t1411 - t1406 + t1445) * t1303) * t1339) * t1526 + t1303 * t1341 * t1481 - ((t1339 * t1395 * t1525 + t1321 * t1373) * t1389 + qJ(2,3) * t1321 * t1490) * t1374 * t1517 - t1315 * t1339 * t1386 + (-t1289 * t1388 - t1292 * t1387) * t1339 * t1373;
t1444 = -t1362 - t1448;
t1300 = t1303 ^ 2;
t1532 = t1300 * t1340;
t1436 = qJ(2,3) * t1532 - t1444;
t1567 = -t1243 * t1404 + t1327 - t1436;
t1363 = g(3) * t1392;
t1338 = qJ(2,2) * t1398 - t1372 * t1392;
t1502 = t1357 * t1397;
t1290 = (pkin(3) * t1502 - t1338 * t1360) * t1391 + t1360 * t1465 + qJ(2,2) * t1502;
t1496 = t1360 * t1397;
t1293 = (pkin(3) * t1496 + t1338 * t1357) * t1391 + t1563 * t1357 * t1543 + qJ(2,2) * t1496;
t1316 = t1350 * t1392 + t1372 * t1398;
t1407 = qJ(2,2) ^ 2;
t1492 = t1372 * t1375;
t1480 = t1259 * t1492;
t1489 = t1375 * t1412;
t1516 = t1323 * t1342;
t1447 = (-t1486 * t1492 + (t1480 + (-0.2e1 * qJ(2,2) * t1546 + t1380 * t1411 - t1407 + t1445) * t1304) * t1342) * t1524 + t1304 * t1344 * t1480 - ((t1342 * t1397 * t1523 + t1323 * t1375) * t1391 + qJ(2,2) * t1323 * t1489) * t1376 * t1516 - t1316 * t1342 * t1386 + (-t1290 * t1388 - t1293 * t1387) * t1342 * t1375;
t1443 = -t1363 - t1447;
t1301 = t1304 ^ 2;
t1530 = t1301 * t1343;
t1437 = qJ(2,2) * t1530 - t1443;
t1566 = -t1244 * t1404 + t1328 - t1437;
t1364 = g(3) * t1394;
t1336 = qJ(2,1) * t1400 - t1372 * t1394;
t1500 = t1358 * t1399;
t1288 = (pkin(3) * t1500 - t1336 * t1361) * t1393 + t1361 * t1464 + qJ(2,1) * t1500;
t1494 = t1361 * t1399;
t1291 = (pkin(3) * t1494 + t1336 * t1358) * t1393 + t1562 * t1358 * t1542 + qJ(2,1) * t1494;
t1317 = t1351 * t1394 + t1372 * t1400;
t1408 = qJ(2,1) ^ 2;
t1491 = t1372 * t1377;
t1479 = t1260 * t1491;
t1488 = t1377 * t1412;
t1515 = t1325 * t1345;
t1446 = (-t1485 * t1491 + (t1479 + (-0.2e1 * qJ(2,1) * t1545 + t1381 * t1411 - t1408 + t1445) * t1305) * t1345) * t1522 + t1305 * t1347 * t1479 - ((t1345 * t1399 * t1521 + t1325 * t1377) * t1393 + qJ(2,1) * t1325 * t1488) * t1378 * t1515 - t1317 * t1345 * t1386 + (-t1288 * t1388 - t1291 * t1387) * t1345 * t1377;
t1442 = -t1364 - t1446;
t1302 = t1305 ^ 2;
t1528 = t1302 * t1346;
t1438 = qJ(2,1) * t1528 - t1442;
t1565 = -t1245 * t1404 + t1329 - t1438;
t1561 = t1472 * t1395;
t1560 = t1470 * t1397;
t1559 = t1468 * t1399;
t1552 = 0.2e1 * qJ(2,3);
t1231 = t1243 * t1552 + t1573;
t1228 = t1404 * t1472 + t1231;
t1297 = -t1373 * t1561 + (-t1356 * t1388 - t1359 * t1387) * t1490;
t1460 = t1303 * t1412 * t1517;
t1457 = t1373 * t1460;
t1272 = t1297 * t1404 + t1457 * t1552;
t1279 = t1297 * t1395 - t1373 * t1520;
t1353 = 0.2e1 * t1379 - 0.1e1;
t1405 = pkin(1) * g(3);
t1535 = t1297 * t1389;
t1538 = t1243 * t1395;
t1555 = 2 * MDP(8);
t1429 = MDP(7) * (0.2e1 * t1460 + t1538) * t1395 + MDP(1) * t1243 + MDP(10) * (-t1535 - t1561) + MDP(12) * (t1228 * t1389 - t1272 * t1395) + MDP(13) * (t1228 * t1395 + t1272 * t1389) + MDP(2) * (t1362 - t1327) + MDP(3) * t1312 + MDP(4) * (t1327 + t1444 - 0.2e1 * t1550) + MDP(5) * t1231 + MDP(6) * (t1243 * t1406 + t1405 * t1390 + t1573 * qJ(2,3) + (t1448 - t1570) * pkin(1)) + MDP(9) * t1279 + (t1353 * t1457 - t1389 * t1538) * t1555;
t1558 = t1390 * t1429;
t1553 = 0.2e1 * qJ(2,2);
t1232 = t1244 * t1553 + t1572;
t1229 = t1404 * t1470 + t1232;
t1298 = -t1375 * t1560 + (-t1357 * t1388 - t1360 * t1387) * t1489;
t1459 = t1304 * t1412 * t1516;
t1456 = t1375 * t1459;
t1270 = t1298 * t1404 + t1456 * t1553;
t1280 = t1298 * t1397 - t1375 * t1519;
t1354 = 0.2e1 * t1380 - 0.1e1;
t1534 = t1298 * t1391;
t1537 = t1244 * t1397;
t1428 = MDP(7) * (0.2e1 * t1459 + t1537) * t1397 + MDP(1) * t1244 + MDP(10) * (-t1534 - t1560) + MDP(12) * (t1229 * t1391 - t1270 * t1397) + MDP(13) * (t1229 * t1397 + t1270 * t1391) + MDP(2) * (t1363 - t1328) + MDP(3) * t1313 + MDP(4) * (t1328 + t1443 - 0.2e1 * t1549) + MDP(5) * t1232 + MDP(6) * (t1244 * t1407 + t1405 * t1392 + t1572 * qJ(2,2) + (t1447 - t1569) * pkin(1)) + MDP(9) * t1280 + (t1354 * t1456 - t1391 * t1537) * t1555;
t1557 = t1392 * t1428;
t1554 = 0.2e1 * qJ(2,1);
t1233 = t1245 * t1554 + t1571;
t1230 = t1404 * t1468 + t1233;
t1299 = -t1377 * t1559 + (-t1358 * t1388 - t1361 * t1387) * t1488;
t1458 = t1305 * t1412 * t1515;
t1455 = t1377 * t1458;
t1271 = t1299 * t1404 + t1455 * t1554;
t1281 = t1299 * t1399 - t1377 * t1518;
t1355 = 0.2e1 * t1381 - 0.1e1;
t1533 = t1299 * t1393;
t1536 = t1245 * t1399;
t1427 = MDP(7) * (0.2e1 * t1458 + t1536) * t1399 + MDP(1) * t1245 + MDP(10) * (-t1533 - t1559) + MDP(12) * (t1230 * t1393 - t1271 * t1399) + MDP(13) * (t1230 * t1399 + t1271 * t1393) + MDP(2) * (t1364 - t1329) + MDP(3) * t1314 + MDP(4) * (t1329 + t1442 - 0.2e1 * t1548) + MDP(5) * t1233 + MDP(6) * (t1245 * t1408 + t1405 * t1394 + t1571 * qJ(2,1) + (t1446 - t1568) * pkin(1)) + MDP(9) * t1281 + (t1355 * t1455 - t1393 * t1536) * t1555;
t1556 = t1394 * t1427;
t1531 = t1300 * t1341;
t1529 = t1301 * t1344;
t1527 = t1302 * t1347;
t1505 = t1356 * t1373;
t1503 = t1357 * t1375;
t1501 = t1358 * t1377;
t1499 = t1359 * t1373;
t1497 = t1360 * t1375;
t1495 = t1361 * t1377;
t1484 = t1373 * t1538;
t1483 = t1375 * t1537;
t1482 = t1377 * t1536;
t1478 = t1395 * t1532;
t1477 = t1373 * t1531;
t1476 = t1397 * t1530;
t1475 = t1375 * t1529;
t1474 = t1399 * t1528;
t1473 = t1377 * t1527;
t1463 = t1300 * t1353 * t1510;
t1462 = t1301 * t1354 * t1508;
t1461 = t1302 * t1355 * t1506;
t1441 = MDP(12) * (-t1389 * t1532 + t1279) + MDP(13) * (-t1535 + (-t1472 - t1532) * t1395) + MDP(4) * t1243 + MDP(6) * (-t1436 + t1570);
t1440 = MDP(12) * (-t1391 * t1530 + t1280) + MDP(13) * (-t1534 + (-t1470 - t1530) * t1397) + MDP(4) * t1244 + MDP(6) * (-t1437 + t1569);
t1439 = MDP(12) * (-t1393 * t1528 + t1281) + MDP(13) * (-t1533 + (-t1468 - t1528) * t1399) + MDP(4) * t1245 + MDP(6) * (-t1438 + t1568);
t1435 = t1441 * t1373;
t1434 = t1440 * t1375;
t1433 = t1439 * t1377;
t1332 = g(1) * t1358 + g(2) * t1361;
t1331 = g(1) * t1357 + g(2) * t1360;
t1330 = g(1) * t1356 + g(2) * t1359;
t1215 = t1332 * t1399 - t1565 * t1393;
t1214 = t1332 * t1393 + t1565 * t1399;
t1213 = t1331 * t1397 - t1566 * t1391;
t1212 = t1331 * t1391 + t1566 * t1397;
t1211 = t1330 * t1395 - t1567 * t1389;
t1210 = t1330 * t1389 + t1567 * t1395;
t1 = [(-t1288 * t1473 - t1289 * t1477 - t1290 * t1475) * MDP(5) + (t1388 - g(1)) * MDP(14) + (t1288 * t1433 + t1361 * t1556) * t1345 + (t1290 * t1434 + t1360 * t1557) * t1342 + (t1289 * t1435 + t1359 * t1558) * t1339 + ((-t1356 * t1478 - t1357 * t1476 - t1358 * t1474) * MDP(7) + (-t1356 * t1463 - t1357 * t1462 - t1358 * t1461) * MDP(8) + (-t1356 * t1484 - t1357 * t1483 - t1358 * t1482) * MDP(9) + (t1243 * t1356 + t1244 * t1357 + t1245 * t1358) * MDP(10) + (-t1297 * t1505 - t1298 * t1503 - t1299 * t1501) * MDP(11) + (-t1210 * t1505 - t1212 * t1503 - t1214 * t1501) * MDP(12) + (-t1211 * t1505 - t1213 * t1503 - t1215 * t1501) * MDP(13)) * t1412; (-t1291 * t1473 - t1292 * t1477 - t1293 * t1475) * MDP(5) + (t1387 - g(2)) * MDP(14) + (t1291 * t1433 - t1358 * t1556) * t1345 + (t1293 * t1434 - t1357 * t1557) * t1342 + (t1292 * t1435 - t1356 * t1558) * t1339 + ((-t1359 * t1478 - t1360 * t1476 - t1361 * t1474) * MDP(7) + (-t1359 * t1463 - t1360 * t1462 - t1361 * t1461) * MDP(8) + (-t1359 * t1484 - t1360 * t1483 - t1361 * t1482) * MDP(9) + (t1243 * t1359 + t1244 * t1360 + t1245 * t1361) * MDP(10) + (-t1297 * t1499 - t1298 * t1497 - t1299 * t1495) * MDP(11) + (-t1210 * t1499 - t1212 * t1497 - t1214 * t1495) * MDP(12) + (-t1211 * t1499 - t1213 * t1497 - t1215 * t1495) * MDP(13)) * t1412; (-t1315 * t1531 - t1316 * t1529 - t1317 * t1527) * MDP(5) + (t1386 - g(3)) * MDP(14) + (t1439 * t1317 + t1427 * t1400) * t1345 + (t1440 * t1316 + t1428 * t1398) * t1342 + (t1441 * t1315 + t1429 * t1396) * t1339;];
tauX  = t1;
