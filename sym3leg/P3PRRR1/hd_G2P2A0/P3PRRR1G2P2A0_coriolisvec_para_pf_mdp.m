% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR1G2P2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G2P2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G2P2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:35
% EndTime: 2020-03-09 21:18:37
% DurationCPUTime: 2.55s
% Computational Cost: add. (19377->168), mult. (11664->366), div. (2469->9), fcn. (11733->24), ass. (0->177)
t1427 = pkin(7) + qJ(2,3);
t1418 = qJ(3,3) + t1427;
t1406 = sin(t1418);
t1409 = cos(t1418);
t1412 = sin(t1427);
t1415 = cos(t1427);
t1388 = -t1415 * t1406 + t1409 * t1412;
t1380 = 0.1e1 / t1388 ^ 2;
t1430 = legFrame(3,2);
t1421 = sin(t1430);
t1424 = cos(t1430);
t1440 = xDP(2);
t1441 = xDP(1);
t1394 = -t1421 * t1440 + t1424 * t1441;
t1439 = xDP(3);
t1460 = -t1394 * t1406 - t1409 * t1439;
t1554 = t1460 ^ 2 * t1380;
t1428 = pkin(7) + qJ(2,2);
t1419 = qJ(3,2) + t1428;
t1407 = sin(t1419);
t1410 = cos(t1419);
t1413 = sin(t1428);
t1416 = cos(t1428);
t1389 = -t1416 * t1407 + t1410 * t1413;
t1383 = 0.1e1 / t1389 ^ 2;
t1431 = legFrame(2,2);
t1422 = sin(t1431);
t1425 = cos(t1431);
t1395 = -t1422 * t1440 + t1425 * t1441;
t1459 = -t1395 * t1407 - t1410 * t1439;
t1553 = t1459 ^ 2 * t1383;
t1429 = pkin(7) + qJ(2,1);
t1420 = qJ(3,1) + t1429;
t1408 = sin(t1420);
t1411 = cos(t1420);
t1414 = sin(t1429);
t1417 = cos(t1429);
t1390 = -t1417 * t1408 + t1411 * t1414;
t1386 = 0.1e1 / t1390 ^ 2;
t1432 = legFrame(1,2);
t1423 = sin(t1432);
t1426 = cos(t1432);
t1396 = -t1423 * t1440 + t1426 * t1441;
t1458 = -t1396 * t1408 - t1411 * t1439;
t1552 = t1458 ^ 2 * t1386;
t1551 = 0.1e1 / t1388;
t1550 = 0.1e1 / t1389;
t1549 = 0.1e1 / t1390;
t1444 = 0.1e1 / pkin(2);
t1525 = t1551 * t1444;
t1364 = -pkin(2) * (t1394 * t1412 + t1415 * t1439) + t1460 * pkin(3);
t1443 = 0.1e1 / pkin(3);
t1539 = t1364 * t1443;
t1359 = (-t1460 + t1539) * t1525;
t1442 = pkin(3) ^ 2;
t1457 = t1406 * t1412 + t1409 * t1415;
t1512 = 0.2e1 * pkin(3) * t1444;
t1533 = t1460 * t1551;
t1548 = (t1359 * t1442 + (-t1533 + t1457 * (-t1460 + t1539 / 0.2e1) * t1551 * t1512) * pkin(2)) * t1443;
t1522 = t1550 * t1444;
t1365 = -pkin(2) * (t1395 * t1413 + t1416 * t1439) + t1459 * pkin(3);
t1538 = t1365 * t1443;
t1361 = (-t1459 + t1538) * t1522;
t1456 = t1407 * t1413 + t1410 * t1416;
t1532 = t1459 * t1550;
t1547 = (t1361 * t1442 + (-t1532 + t1456 * (-t1459 + t1538 / 0.2e1) * t1550 * t1512) * pkin(2)) * t1443;
t1519 = t1549 * t1444;
t1366 = -pkin(2) * (t1396 * t1414 + t1417 * t1439) + t1458 * pkin(3);
t1537 = t1366 * t1443;
t1363 = (-t1458 + t1537) * t1519;
t1455 = t1408 * t1414 + t1411 * t1417;
t1531 = t1458 * t1549;
t1546 = (t1363 * t1442 + (-t1531 + t1455 * (-t1458 + t1537 / 0.2e1) * t1549 * t1512) * pkin(2)) * t1443;
t1545 = (-0.2e1 * t1460 + t1539) * t1525 * t1364;
t1544 = t1359 * t1364;
t1543 = (-0.2e1 * t1459 + t1538) * t1522 * t1365;
t1542 = t1361 * t1365;
t1541 = (-0.2e1 * t1458 + t1537) * t1519 * t1366;
t1540 = t1363 * t1366;
t1536 = t1551 * t1554;
t1535 = t1550 * t1553;
t1534 = t1549 * t1552;
t1530 = t1551 ^ 2;
t1529 = t1550 ^ 2;
t1528 = t1549 ^ 2;
t1400 = -pkin(2) * t1415 - pkin(3) * t1409;
t1527 = t1551 * t1400;
t1526 = t1551 * t1409;
t1401 = -pkin(2) * t1416 - pkin(3) * t1410;
t1524 = t1550 * t1401;
t1523 = t1550 * t1410;
t1402 = -pkin(2) * t1417 - pkin(3) * t1411;
t1521 = t1549 * t1402;
t1520 = t1549 * t1411;
t1397 = pkin(2) * t1412 + pkin(3) * t1406;
t1518 = t1397 * t1421;
t1398 = pkin(2) * t1413 + pkin(3) * t1407;
t1517 = t1398 * t1422;
t1399 = pkin(2) * t1414 + pkin(3) * t1408;
t1516 = t1399 * t1423;
t1515 = t1406 * t1424;
t1514 = t1407 * t1425;
t1513 = t1408 * t1426;
t1352 = -pkin(3) * t1359 + t1457 * t1533;
t1445 = 0.1e1 / pkin(2) ^ 2;
t1346 = (t1352 * t1460 + t1544) * t1445 * t1530;
t1511 = t1346 * t1527;
t1353 = -pkin(3) * t1361 + t1456 * t1532;
t1347 = (t1353 * t1459 + t1542) * t1445 * t1529;
t1510 = t1347 * t1524;
t1354 = -pkin(3) * t1363 + t1455 * t1531;
t1348 = (t1354 * t1458 + t1540) * t1445 * t1528;
t1509 = t1348 * t1521;
t1508 = t1380 * t1545;
t1507 = t1383 * t1543;
t1506 = t1386 * t1541;
t1505 = t1400 * t1536;
t1504 = t1401 * t1535;
t1503 = t1402 * t1534;
t1502 = t1551 * t1397 * t1424;
t1501 = t1551 * t1406 * t1421;
t1500 = t1550 * t1398 * t1425;
t1499 = t1550 * t1407 * t1422;
t1498 = t1549 * t1399 * t1426;
t1497 = t1549 * t1408 * t1423;
t1496 = t1551 * t1518;
t1495 = t1551 * t1515;
t1494 = t1550 * t1517;
t1493 = t1550 * t1514;
t1492 = t1549 * t1516;
t1491 = t1549 * t1513;
t1468 = t1359 * (t1457 * pkin(2) + pkin(3)) * t1380 * t1539;
t1341 = (t1468 - (0.2e1 * t1544 - (-0.2e1 * t1352 - t1548) * t1460) * t1530) * t1445;
t1490 = t1341 * t1526;
t1466 = t1361 * (t1456 * pkin(2) + pkin(3)) * t1383 * t1538;
t1343 = (t1466 - (0.2e1 * t1542 - (-0.2e1 * t1353 - t1547) * t1459) * t1529) * t1445;
t1489 = t1343 * t1523;
t1464 = t1363 * (t1455 * pkin(2) + pkin(3)) * t1386 * t1537;
t1345 = (t1464 - (0.2e1 * t1540 - (-0.2e1 * t1354 - t1546) * t1458) * t1528) * t1445;
t1488 = t1345 * t1520;
t1487 = t1409 * t1508;
t1486 = t1410 * t1507;
t1485 = t1411 * t1506;
t1484 = t1518 * t1536;
t1483 = t1517 * t1535;
t1482 = t1516 * t1534;
t1481 = t1346 * t1496;
t1480 = t1341 * t1501;
t1479 = t1347 * t1494;
t1478 = t1343 * t1499;
t1477 = t1348 * t1492;
t1476 = t1345 * t1497;
t1475 = t1346 * t1502;
t1474 = t1341 * t1495;
t1473 = t1347 * t1500;
t1472 = t1343 * t1493;
t1471 = t1348 * t1498;
t1470 = t1345 * t1491;
t1469 = t1508 * t1515;
t1467 = t1507 * t1514;
t1465 = t1506 * t1513;
t1463 = t1502 * t1554;
t1462 = t1500 * t1553;
t1461 = t1498 * t1552;
t1454 = t1551 * t1501 * t1545;
t1453 = t1550 * t1499 * t1543;
t1452 = t1549 * t1497 * t1541;
t1438 = cos(qJ(3,1));
t1437 = cos(qJ(3,2));
t1436 = cos(qJ(3,3));
t1435 = sin(qJ(3,1));
t1434 = sin(qJ(3,2));
t1433 = sin(qJ(3,3));
t1344 = (t1464 - (t1540 - (-t1354 - t1546) * t1458) * t1528) * t1445;
t1342 = (t1466 - (t1542 - (-t1353 - t1547) * t1459) * t1529) * t1445;
t1340 = (t1468 - (t1544 - (-t1352 - t1548) * t1460) * t1530) * t1445;
t1 = [(-t1436 * t1474 - t1437 * t1472 - t1438 * t1470) * MDP(6) + (t1433 * t1474 + t1434 * t1472 + t1435 * t1470) * MDP(7) + ((t1346 * t1495 + t1347 * t1493 + t1348 * t1491) * MDP(2) + (-t1340 * t1495 - t1342 * t1493 - t1344 * t1491) * MDP(5)) * t1444 + ((-t1436 * t1475 - t1437 * t1473 - t1438 * t1471) * MDP(6) + (t1433 * t1475 + t1434 * t1473 + t1435 * t1471) * MDP(7) + ((t1433 * t1463 + t1434 * t1462 + t1435 * t1461) * MDP(6) + (t1436 * t1463 + t1437 * t1462 + t1438 * t1461) * MDP(7)) * t1445 + ((t1340 * t1502 + t1342 * t1500 + t1344 * t1498) * MDP(5) + (t1433 * t1469 + t1434 * t1467 + t1435 * t1465) * MDP(6) + (t1436 * t1469 + t1437 * t1467 + t1438 * t1465) * MDP(7)) * t1444) * t1443; (t1436 * t1480 + t1437 * t1478 + t1438 * t1476) * MDP(6) + (-t1433 * t1480 - t1434 * t1478 - t1435 * t1476) * MDP(7) + ((-t1346 * t1501 - t1347 * t1499 - t1348 * t1497) * MDP(2) + (t1340 * t1501 + t1342 * t1499 + t1344 * t1497) * MDP(5)) * t1444 + ((t1436 * t1481 + t1437 * t1479 + t1438 * t1477) * MDP(6) + (-t1433 * t1481 - t1434 * t1479 - t1435 * t1477) * MDP(7) + ((-t1433 * t1484 - t1434 * t1483 - t1435 * t1482) * MDP(6) + (-t1436 * t1484 - t1437 * t1483 - t1438 * t1482) * MDP(7)) * t1445 + ((-t1340 * t1496 - t1342 * t1494 - t1344 * t1492) * MDP(5) + (-t1433 * t1454 - t1434 * t1453 - t1435 * t1452) * MDP(6) + (-t1436 * t1454 - t1437 * t1453 - t1438 * t1452) * MDP(7)) * t1444) * t1443; (-t1436 * t1490 - t1437 * t1489 - t1438 * t1488) * MDP(6) + (t1433 * t1490 + t1434 * t1489 + t1435 * t1488) * MDP(7) + ((t1346 * t1526 + t1347 * t1523 + t1348 * t1520) * MDP(2) + (-t1340 * t1526 - t1342 * t1523 - t1344 * t1520) * MDP(5)) * t1444 + ((t1436 * t1511 + t1437 * t1510 + t1438 * t1509) * MDP(6) + (-t1433 * t1511 - t1434 * t1510 - t1435 * t1509) * MDP(7) + ((-t1433 * t1505 - t1434 * t1504 - t1435 * t1503) * MDP(6) + (-t1436 * t1505 - t1437 * t1504 - t1438 * t1503) * MDP(7)) * t1445 + ((-t1340 * t1527 - t1342 * t1524 - t1344 * t1521) * MDP(5) + (t1433 * t1487 + t1434 * t1486 + t1435 * t1485) * MDP(6) + (t1436 * t1487 + t1437 * t1486 + t1438 * t1485) * MDP(7)) * t1444) * t1443;];
taucX  = t1;
