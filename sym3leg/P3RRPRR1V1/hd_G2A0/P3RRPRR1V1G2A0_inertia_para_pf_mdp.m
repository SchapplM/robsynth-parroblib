% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR1V1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RRPRR1V1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:34:13
% EndTime: 2020-08-06 19:34:15
% DurationCPUTime: 2.21s
% Computational Cost: add. (1569->264), mult. (2315->555), div. (855->14), fcn. (2499->18), ass. (0->238)
t1503 = 2 * MDP(11);
t1401 = legFrame(3,2);
t1372 = sin(t1401);
t1375 = cos(t1401);
t1404 = sin(qJ(2,3));
t1405 = sin(qJ(1,3));
t1410 = cos(qJ(2,3));
t1513 = t1405 * t1410;
t1363 = -t1372 * t1513 + t1404 * t1375;
t1366 = t1404 * t1372 + t1375 * t1513;
t1561 = t1363 * t1366;
t1402 = legFrame(2,2);
t1373 = sin(t1402);
t1376 = cos(t1402);
t1406 = sin(qJ(2,2));
t1407 = sin(qJ(1,2));
t1412 = cos(qJ(2,2));
t1510 = t1407 * t1412;
t1364 = -t1373 * t1510 + t1406 * t1376;
t1367 = t1406 * t1373 + t1376 * t1510;
t1560 = t1364 * t1367;
t1403 = legFrame(1,2);
t1374 = sin(t1403);
t1377 = cos(t1403);
t1408 = sin(qJ(2,1));
t1409 = sin(qJ(1,1));
t1414 = cos(qJ(2,1));
t1507 = t1409 * t1414;
t1365 = -t1374 * t1507 + t1408 * t1377;
t1368 = t1408 * t1374 + t1377 * t1507;
t1559 = t1365 * t1368;
t1558 = 2 * MDP(5);
t1557 = 2 * MDP(6);
t1556 = 2 * MDP(7);
t1416 = pkin(1) + pkin(2);
t1396 = 1 / t1416;
t1555 = pkin(1) * t1396;
t1554 = pkin(1) * t1410;
t1553 = pkin(1) * t1412;
t1552 = pkin(1) * t1414;
t1551 = MDP(6) * t1396;
t1550 = MDP(7) * t1396;
t1549 = MDP(8) / t1416 ^ 2;
t1398 = pkin(3) + qJ(3,3);
t1411 = cos(qJ(1,3));
t1506 = t1410 * t1416;
t1433 = -t1398 * t1411 + t1405 * t1506;
t1514 = t1404 * t1416;
t1348 = -t1433 * t1372 + t1375 * t1514;
t1378 = 0.1e1 / t1398;
t1548 = t1348 * t1378;
t1400 = pkin(3) + qJ(3,1);
t1415 = cos(qJ(1,1));
t1504 = t1414 * t1416;
t1431 = -t1400 * t1415 + t1409 * t1504;
t1508 = t1408 * t1416;
t1349 = -t1431 * t1374 + t1377 * t1508;
t1382 = 0.1e1 / t1400;
t1547 = t1349 * t1382;
t1399 = pkin(3) + qJ(3,2);
t1413 = cos(qJ(1,2));
t1505 = t1412 * t1416;
t1432 = -t1399 * t1413 + t1407 * t1505;
t1511 = t1406 * t1416;
t1350 = -t1432 * t1373 + t1376 * t1511;
t1380 = 0.1e1 / t1399;
t1546 = t1350 * t1380;
t1351 = t1372 * t1514 + t1433 * t1375;
t1545 = t1351 * t1378;
t1352 = t1373 * t1511 + t1432 * t1376;
t1544 = t1352 * t1380;
t1353 = t1374 * t1508 + t1431 * t1377;
t1543 = t1353 * t1382;
t1369 = t1405 * t1398 + t1411 * t1506;
t1354 = (-t1411 * t1554 + t1369) * t1378;
t1542 = t1354 * t1378;
t1370 = t1407 * t1399 + t1413 * t1505;
t1355 = (-t1413 * t1553 + t1370) * t1380;
t1541 = t1355 * t1380;
t1371 = t1409 * t1400 + t1415 * t1504;
t1356 = (-t1415 * t1552 + t1371) * t1382;
t1540 = t1356 * t1382;
t1387 = 0.1e1 / t1410;
t1539 = t1372 * t1387;
t1390 = 0.1e1 / t1412;
t1538 = t1373 * t1390;
t1393 = 0.1e1 / t1414;
t1537 = t1374 * t1393;
t1536 = t1375 * t1387;
t1535 = t1376 * t1390;
t1534 = t1377 * t1393;
t1533 = t1378 * t1387;
t1532 = t1378 * t1404;
t1531 = t1378 * t1411;
t1379 = 0.1e1 / t1398 ^ 2;
t1530 = t1379 * t1387;
t1421 = t1410 ^ 2;
t1388 = 0.1e1 / t1421;
t1529 = t1379 * t1388;
t1528 = t1380 * t1390;
t1527 = t1380 * t1406;
t1526 = t1380 * t1413;
t1381 = 0.1e1 / t1399 ^ 2;
t1525 = t1381 * t1390;
t1422 = t1412 ^ 2;
t1391 = 0.1e1 / t1422;
t1524 = t1381 * t1391;
t1523 = t1382 * t1393;
t1522 = t1382 * t1408;
t1521 = t1382 * t1415;
t1383 = 0.1e1 / t1400 ^ 2;
t1520 = t1383 * t1393;
t1423 = t1414 ^ 2;
t1394 = 0.1e1 / t1423;
t1519 = t1383 * t1394;
t1518 = t1411 ^ 2 * t1379;
t1517 = t1413 ^ 2 * t1381;
t1516 = t1415 ^ 2 * t1383;
t1515 = t1404 * t1411;
t1512 = t1406 * t1413;
t1509 = t1408 * t1415;
t1502 = 0.2e1 * qJ(3,1) * t1382;
t1501 = 0.2e1 * qJ(3,2) * t1380;
t1500 = 0.2e1 * qJ(3,3) * t1378;
t1499 = t1404 * t1555;
t1498 = t1406 * t1555;
t1497 = t1408 * t1555;
t1496 = qJ(3,1) * t1522;
t1495 = qJ(3,2) * t1527;
t1494 = qJ(3,3) * t1532;
t1493 = t1363 * t1533;
t1492 = t1364 * t1528;
t1491 = t1365 * t1523;
t1490 = t1366 * t1533;
t1489 = t1367 * t1528;
t1488 = t1368 * t1523;
t1417 = qJ(3,3) ^ 2;
t1420 = pkin(1) ^ 2;
t1333 = (-t1369 * t1554 + (t1420 * t1421 + t1417) * t1411) * t1378;
t1487 = t1333 * t1533;
t1486 = t1388 * t1532;
t1384 = t1404 ^ 2;
t1485 = t1384 * t1529;
t1484 = t1404 * t1530;
t1483 = t1411 * t1530;
t1482 = t1379 * t1515;
t1418 = qJ(3,2) ^ 2;
t1334 = (-t1370 * t1553 + (t1420 * t1422 + t1418) * t1413) * t1380;
t1481 = t1334 * t1528;
t1480 = t1391 * t1527;
t1385 = t1406 ^ 2;
t1479 = t1385 * t1524;
t1478 = t1406 * t1525;
t1477 = t1413 * t1525;
t1476 = t1381 * t1512;
t1419 = qJ(3,1) ^ 2;
t1335 = (-t1371 * t1552 + (t1420 * t1423 + t1419) * t1415) * t1382;
t1475 = t1335 * t1523;
t1474 = t1394 * t1522;
t1386 = t1408 ^ 2;
t1473 = t1386 * t1519;
t1472 = t1408 * t1520;
t1471 = t1415 * t1520;
t1470 = t1383 * t1509;
t1469 = t1363 * t1483;
t1468 = t1364 * t1477;
t1467 = t1365 * t1471;
t1466 = t1366 * t1483;
t1465 = t1529 * t1561;
t1464 = t1367 * t1477;
t1463 = t1524 * t1560;
t1462 = t1368 * t1471;
t1461 = t1519 * t1559;
t1460 = t1372 * t1486;
t1459 = t1373 * t1480;
t1458 = t1374 * t1474;
t1457 = t1375 * t1486;
t1456 = t1376 * t1480;
t1455 = t1377 * t1474;
t1454 = t1515 * t1533;
t1453 = t1384 * t1483;
t1452 = t1512 * t1528;
t1451 = t1385 * t1477;
t1450 = t1509 * t1523;
t1449 = t1386 * t1471;
t1448 = qJ(3,1) * t1393 * t1497;
t1447 = qJ(3,2) * t1390 * t1498;
t1446 = qJ(3,3) * t1387 * t1499;
t1445 = t1372 * t1454;
t1444 = t1373 * t1452;
t1443 = t1374 * t1450;
t1442 = t1375 * t1454;
t1441 = t1376 * t1452;
t1440 = t1377 * t1450;
t1424 = t1440 + t1441 + t1442;
t1439 = (t1363 * t1453 + t1364 * t1451 + t1365 * t1449) * MDP(4) + (t1363 * t1482 + t1364 * t1476 + t1365 * t1470) * t1558 + (t1467 + t1468 + t1469) * MDP(1) + t1424 * t1551 + (t1375 * t1531 + t1376 * t1526 + t1377 * t1521) * t1550;
t1425 = t1443 + t1444 + t1445;
t1438 = (t1366 * t1453 + t1367 * t1451 + t1368 * t1449) * MDP(4) + (t1366 * t1482 + t1367 * t1476 + t1368 * t1470) * t1558 + (t1462 + t1464 + t1466) * MDP(1) + t1425 * t1551 + (t1372 * t1531 + t1373 * t1526 + t1374 * t1521) * t1550;
t1437 = t1387 * t1417 + t1410 * t1420;
t1436 = t1390 * t1418 + t1412 * t1420;
t1435 = t1393 * t1419 + t1414 * t1420;
t1428 = t1382 * (t1365 * t1374 + t1368 * t1377);
t1429 = t1380 * (t1364 * t1373 + t1367 * t1376);
t1430 = t1378 * (t1363 * t1372 + t1366 * t1375);
t1434 = (t1404 * t1388 * t1430 + t1406 * t1391 * t1429 + t1408 * t1394 * t1428) * t1551 + (t1387 * t1430 + t1390 * t1429 + t1393 * t1428) * t1550 + (t1472 * t1559 + t1478 * t1560 + t1484 * t1561) * t1558 + (t1384 * t1465 + t1385 * t1463 + t1386 * t1461) * MDP(4) + (t1461 + t1463 + t1465) * MDP(1) + (t1372 * t1375 * t1388 + t1373 * t1376 * t1391 + t1374 * t1377 * t1394) * t1549;
t1427 = t1363 * t1457 + t1364 * t1456 + t1365 * t1455;
t1426 = t1366 * t1460 + t1367 * t1459 + t1368 * t1458;
t1362 = t1368 ^ 2;
t1361 = t1367 ^ 2;
t1360 = t1366 ^ 2;
t1359 = t1365 ^ 2;
t1358 = t1364 ^ 2;
t1357 = t1363 ^ 2;
t1347 = (t1368 * t1502 - t1374 * t1497) * t1393;
t1346 = (t1367 * t1501 - t1373 * t1498) * t1390;
t1345 = (t1366 * t1500 - t1372 * t1499) * t1387;
t1344 = (t1365 * t1502 - t1377 * t1497) * t1393;
t1343 = (t1364 * t1501 - t1376 * t1498) * t1390;
t1342 = (t1363 * t1500 - t1375 * t1499) * t1387;
t1341 = (-t1368 * t1496 + t1374 * t1555) * t1393;
t1340 = (-t1365 * t1496 + t1377 * t1555) * t1393;
t1339 = (-t1367 * t1495 + t1373 * t1555) * t1390;
t1338 = (-t1364 * t1495 + t1376 * t1555) * t1390;
t1337 = (-t1366 * t1494 + t1372 * t1555) * t1387;
t1336 = (-t1363 * t1494 + t1375 * t1555) * t1387;
t1327 = (-pkin(1) * t1368 + t1353) * t1382;
t1326 = (-pkin(1) * t1367 + t1352) * t1380;
t1325 = (-pkin(1) * t1366 + t1351) * t1378;
t1324 = (-pkin(1) * t1365 + t1349) * t1382;
t1323 = (-pkin(1) * t1364 + t1350) * t1380;
t1322 = (-pkin(1) * t1363 + t1348) * t1378;
t1315 = -t1374 * t1448 + (-t1353 * t1552 + t1435 * t1368) * t1382;
t1314 = -t1377 * t1448 + (-t1349 * t1552 + t1435 * t1365) * t1382;
t1313 = -t1373 * t1447 + (-t1352 * t1553 + t1436 * t1367) * t1380;
t1312 = -t1376 * t1447 + (-t1350 * t1553 + t1436 * t1364) * t1380;
t1311 = -t1372 * t1446 + (-t1351 * t1554 + t1437 * t1366) * t1378;
t1310 = -t1375 * t1446 + (-t1348 * t1554 + t1437 * t1363) * t1378;
t1 = [(t1360 * t1529 + t1361 * t1524 + t1362 * t1519) * MDP(1) + (t1360 * t1485 + t1361 * t1479 + t1362 * t1473) * MDP(4) + (t1360 * t1484 + t1361 * t1478 + t1362 * t1472) * t1558 + (t1345 * t1490 + t1346 * t1489 + t1347 * t1488) * MDP(11) + (t1311 * t1490 + t1313 * t1489 + t1315 * t1488 + t1325 * t1545 + t1326 * t1544 + t1327 * t1543) * MDP(12) + MDP(13) + (t1372 ^ 2 * t1388 + t1373 ^ 2 * t1391 + t1374 ^ 2 * t1394) * t1549 + (t1426 * t1557 + (t1372 * t1490 + t1373 * t1489 + t1374 * t1488) * t1556 + (-t1426 * MDP(11) + (t1337 * t1539 + t1339 * t1538 + t1341 * t1537) * MDP(12)) * pkin(1)) * t1396; (t1342 * t1490 + t1343 * t1489 + t1344 * t1488) * MDP(11) + (t1310 * t1490 + t1312 * t1489 + t1314 * t1488 + t1322 * t1545 + t1323 * t1544 + t1324 * t1543) * MDP(12) + ((-t1363 * t1460 - t1364 * t1459 - t1365 * t1458) * MDP(11) + (t1336 * t1539 + t1338 * t1538 + t1340 * t1537) * MDP(12)) * t1555 + t1434; (qJ(3,1) * t1462 + qJ(3,2) * t1464 + qJ(3,3) * t1466) * t1503 + (t1351 * t1542 + t1352 * t1541 + t1353 * t1540 + t1366 * t1487 + t1367 * t1481 + t1368 * t1475) * MDP(12) + (-t1425 * MDP(11) + (-qJ(3,1) * t1443 - qJ(3,2) * t1444 - qJ(3,3) * t1445) * MDP(12)) * t1555 + t1438; (t1345 * t1493 + t1346 * t1492 + t1347 * t1491) * MDP(11) + (t1311 * t1493 + t1313 * t1492 + t1315 * t1491 + t1325 * t1548 + t1326 * t1546 + t1327 * t1547) * MDP(12) + ((-t1366 * t1457 - t1367 * t1456 - t1368 * t1455) * MDP(11) + (t1337 * t1536 + t1339 * t1535 + t1341 * t1534) * MDP(12)) * t1555 + t1434; (t1357 * t1529 + t1358 * t1524 + t1359 * t1519) * MDP(1) + (t1357 * t1485 + t1358 * t1479 + t1359 * t1473) * MDP(4) + (t1357 * t1484 + t1358 * t1478 + t1359 * t1472) * t1558 + (t1342 * t1493 + t1343 * t1492 + t1344 * t1491) * MDP(11) + (t1310 * t1493 + t1312 * t1492 + t1314 * t1491 + t1322 * t1548 + t1323 * t1546 + t1324 * t1547) * MDP(12) + MDP(13) + (t1375 ^ 2 * t1388 + t1376 ^ 2 * t1391 + t1377 ^ 2 * t1394) * t1549 + (t1427 * t1557 + (t1375 * t1493 + t1376 * t1492 + t1377 * t1491) * t1556 + (-t1427 * MDP(11) + (t1336 * t1536 + t1338 * t1535 + t1340 * t1534) * MDP(12)) * pkin(1)) * t1396; (qJ(3,1) * t1467 + qJ(3,2) * t1468 + qJ(3,3) * t1469) * t1503 + (t1348 * t1542 + t1349 * t1540 + t1350 * t1541 + t1363 * t1487 + t1364 * t1481 + t1365 * t1475) * MDP(12) + (-t1424 * MDP(11) + (-qJ(3,1) * t1440 - qJ(3,2) * t1441 - qJ(3,3) * t1442) * MDP(12)) * t1555 + t1439; (t1345 * t1531 + t1346 * t1526 + t1347 * t1521) * MDP(11) + ((t1315 * t1415 + t1327 * t1371) * t1382 + (t1313 * t1413 + t1326 * t1370) * t1380 + (t1311 * t1411 + t1325 * t1369) * t1378) * MDP(12) + t1438; (t1342 * t1531 + t1343 * t1526 + t1344 * t1521) * MDP(11) + ((t1314 * t1415 + t1324 * t1371) * t1382 + (t1312 * t1413 + t1323 * t1370) * t1380 + (t1310 * t1411 + t1322 * t1369) * t1378) * MDP(12) + t1439; (t1516 + t1517 + t1518) * MDP(1) + (t1384 * t1518 + t1385 * t1517 + t1386 * t1516) * MDP(4) + ((t1335 * t1415 + t1356 * t1371) * t1382 + (t1334 * t1413 + t1355 * t1370) * t1380 + (t1333 * t1411 + t1354 * t1369) * t1378) * MDP(12) + MDP(13) + (t1404 * t1410 * t1518 + t1406 * t1412 * t1517 + t1408 * t1414 * t1516) * t1558 + (qJ(3,1) * t1516 + qJ(3,2) * t1517 + qJ(3,3) * t1518) * t1503;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;