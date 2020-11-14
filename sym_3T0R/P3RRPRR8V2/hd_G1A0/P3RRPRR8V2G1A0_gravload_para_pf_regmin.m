% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:03:37
% EndTime: 2020-08-06 21:03:38
% DurationCPUTime: 0.83s
% Computational Cost: add. (687->151), mult. (954->243), div. (93->9), fcn. (909->35), ass. (0->115)
t1466 = legFrame(3,3);
t1447 = sin(t1466);
t1450 = cos(t1466);
t1425 = t1447 * g(1) - t1450 * g(2);
t1428 = t1450 * g(1) + t1447 * g(2);
t1473 = sin(qJ(1,3));
t1479 = cos(qJ(1,3));
t1392 = t1473 * t1425 - t1428 * t1479;
t1467 = legFrame(2,3);
t1448 = sin(t1467);
t1451 = cos(t1467);
t1426 = t1448 * g(1) - t1451 * g(2);
t1429 = t1451 * g(1) + t1448 * g(2);
t1475 = sin(qJ(1,2));
t1481 = cos(qJ(1,2));
t1393 = t1475 * t1426 - t1429 * t1481;
t1468 = legFrame(1,3);
t1449 = sin(t1468);
t1452 = cos(t1468);
t1427 = t1449 * g(1) - t1452 * g(2);
t1430 = t1452 * g(1) + t1449 * g(2);
t1477 = sin(qJ(1,1));
t1483 = cos(qJ(1,1));
t1394 = t1477 * t1427 - t1430 * t1483;
t1440 = cos(pkin(7)) * pkin(3) + pkin(2);
t1471 = pkin(5) + qJ(3,1);
t1461 = -pkin(6) - t1471;
t1455 = 0.1e1 / t1461;
t1476 = sin(qJ(2,1));
t1482 = cos(qJ(2,1));
t1512 = sin(pkin(7)) * pkin(3);
t1498 = 0.1e1 / (t1440 * t1482 - t1476 * t1512) * (t1476 * t1440 + t1482 * t1512) * t1455;
t1470 = pkin(5) + qJ(3,2);
t1460 = -pkin(6) - t1470;
t1454 = 0.1e1 / t1460;
t1474 = sin(qJ(2,2));
t1480 = cos(qJ(2,2));
t1499 = 0.1e1 / (t1440 * t1480 - t1474 * t1512) * (t1474 * t1440 + t1480 * t1512) * t1454;
t1469 = pkin(5) + qJ(3,3);
t1459 = -pkin(6) - t1469;
t1453 = 0.1e1 / t1459;
t1472 = sin(qJ(2,3));
t1478 = cos(qJ(2,3));
t1500 = 0.1e1 / (t1440 * t1478 - t1472 * t1512) * (t1472 * t1440 + t1478 * t1512) * t1453;
t1518 = t1392 * t1500 + t1393 * t1499 + t1394 * t1498;
t1424 = t1449 * t1483 + t1452 * t1477;
t1504 = t1424 * t1455;
t1423 = t1448 * t1481 + t1451 * t1475;
t1505 = t1423 * t1454;
t1422 = t1447 * t1479 + t1450 * t1473;
t1506 = t1422 * t1453;
t1517 = t1392 * t1506 + t1393 * t1505 + t1394 * t1504;
t1421 = -t1449 * t1477 + t1452 * t1483;
t1507 = t1421 * t1455;
t1420 = -t1448 * t1475 + t1451 * t1481;
t1508 = t1420 * t1454;
t1419 = -t1447 * t1473 + t1450 * t1479;
t1509 = t1419 * t1453;
t1516 = t1392 * t1509 + t1393 * t1508 + t1394 * t1507;
t1462 = qJ(2,3) + pkin(7);
t1515 = pkin(3) * cos(t1462);
t1463 = qJ(2,2) + pkin(7);
t1514 = pkin(3) * cos(t1463);
t1464 = qJ(2,1) + pkin(7);
t1513 = pkin(3) * cos(t1464);
t1511 = 0.2e1 * pkin(2) * pkin(3);
t1510 = 2 * pkin(1);
t1404 = (t1479 * g(1) + t1473 * g(2)) * t1450 - (t1473 * g(1) - t1479 * g(2)) * t1447;
t1456 = t1478 * pkin(2);
t1431 = 0.1e1 / (t1456 + t1515);
t1503 = t1431 * (-g(3) * t1478 + t1404 * t1472);
t1405 = (t1481 * g(1) + t1475 * g(2)) * t1451 - (t1475 * g(1) - t1481 * g(2)) * t1448;
t1457 = t1480 * pkin(2);
t1432 = 0.1e1 / (t1457 + t1514);
t1502 = t1432 * (-g(3) * t1480 + t1405 * t1474);
t1406 = (g(1) * t1483 + g(2) * t1477) * t1452 - (g(1) * t1477 - g(2) * t1483) * t1449;
t1458 = t1482 * pkin(2);
t1433 = 0.1e1 / (t1458 + t1513);
t1501 = t1433 * (-g(3) * t1482 + t1406 * t1476);
t1395 = t1425 * t1479 + t1428 * t1473;
t1497 = t1395 * t1509;
t1396 = t1426 * t1481 + t1429 * t1475;
t1496 = t1396 * t1508;
t1397 = t1427 * t1483 + t1430 * t1477;
t1495 = t1397 * t1507;
t1494 = t1395 * t1506;
t1493 = t1396 * t1505;
t1492 = t1397 * t1504;
t1491 = t1395 * t1500;
t1490 = t1396 * t1499;
t1489 = t1397 * t1498;
t1488 = pkin(2) ^ 2;
t1487 = pkin(3) ^ 2;
t1486 = 0.2e1 * qJ(2,1);
t1485 = 0.2e1 * qJ(2,2);
t1484 = 0.2e1 * qJ(2,3);
t1443 = t1458 + pkin(1);
t1442 = t1457 + pkin(1);
t1441 = t1456 + pkin(1);
t1439 = t1442 * t1481;
t1438 = t1441 * t1479;
t1437 = t1483 * t1443;
t1436 = t1477 * t1443;
t1435 = t1475 * t1442;
t1434 = t1473 * t1441;
t1418 = -t1475 * t1460 + t1439;
t1417 = -t1473 * t1459 + t1438;
t1416 = -t1477 * t1461 + t1437;
t1415 = t1483 * t1461 + t1436;
t1414 = t1481 * t1460 + t1435;
t1413 = t1479 * t1459 + t1434;
t1388 = t1430 * (-t1471 * t1483 + t1436) + t1427 * (t1471 * t1477 + t1437);
t1387 = t1429 * (-t1470 * t1481 + t1435) + t1426 * (t1470 * t1475 + t1439);
t1386 = t1428 * (-t1469 * t1479 + t1434) + t1425 * (t1469 * t1473 + t1438);
t1 = [0, -t1495 - t1496 - t1497, t1516, 0, 0, 0, 0, 0, -t1478 * t1497 - t1480 * t1496 - t1482 * t1495, t1472 * t1497 + t1474 * t1496 + t1476 * t1495, -t1516, -(t1421 * t1388 - (-t1415 * t1449 + t1416 * t1452 + t1421 * t1513) * t1397) * t1455 - (t1420 * t1387 - (-t1414 * t1448 + t1418 * t1451 + t1420 * t1514) * t1396) * t1454 - (t1419 * t1386 - (-t1413 * t1447 + t1417 * t1450 + t1419 * t1515) * t1395) * t1453, -g(1); 0, -t1492 - t1493 - t1494, t1517, 0, 0, 0, 0, 0, -t1478 * t1494 - t1480 * t1493 - t1482 * t1492, t1472 * t1494 + t1474 * t1493 + t1476 * t1492, -t1517, -(t1424 * t1388 - (t1415 * t1452 + t1416 * t1449 + t1424 * t1513) * t1397) * t1455 - (t1423 * t1387 - (t1414 * t1451 + t1418 * t1448 + t1423 * t1514) * t1396) * t1454 - (t1422 * t1386 - (t1413 * t1450 + t1417 * t1447 + t1422 * t1515) * t1395) * t1453, -g(2); 0, -t1489 - t1490 - t1491, t1518, 0, 0, 0, 0, 0, -t1478 * t1491 - t1480 * t1490 - t1482 * t1489 + t1501 + t1502 + t1503, t1476 * t1489 + t1433 * (g(3) * t1476 + t1406 * t1482) + t1474 * t1490 + t1432 * (g(3) * t1474 + t1405 * t1480) + t1472 * t1491 + t1431 * (g(3) * t1472 + t1404 * t1478), -t1518, -t1388 * t1498 + pkin(2) * t1501 + (sin(t1486 + pkin(7)) * t1511 + t1488 * sin(t1486) + t1487 * sin(0.2e1 * t1464) + (sin(t1464) * pkin(3) + pkin(2) * t1476) * t1510) * t1455 * t1433 * t1397 / 0.2e1 - t1387 * t1499 + pkin(2) * t1502 + (sin(t1485 + pkin(7)) * t1511 + t1488 * sin(t1485) + t1487 * sin(0.2e1 * t1463) + (sin(t1463) * pkin(3) + pkin(2) * t1474) * t1510) * t1454 * t1432 * t1396 / 0.2e1 - t1386 * t1500 + pkin(2) * t1503 + (sin(t1484 + pkin(7)) * t1511 + t1488 * sin(t1484) + t1487 * sin(0.2e1 * t1462) + (sin(t1462) * pkin(3) + pkin(2) * t1472) * t1510) * t1453 * t1431 * t1395 / 0.2e1, -g(3);];
tau_reg  = t1;
