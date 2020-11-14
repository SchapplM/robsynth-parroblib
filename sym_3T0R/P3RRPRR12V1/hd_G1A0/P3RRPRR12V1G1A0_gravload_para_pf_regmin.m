% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V1G1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:02:19
% EndTime: 2020-08-06 19:02:20
% DurationCPUTime: 1.51s
% Computational Cost: add. (909->148), mult. (1635->269), div. (132->6), fcn. (1641->18), ass. (0->148)
t1516 = legFrame(1,3);
t1507 = sin(t1516);
t1510 = cos(t1516);
t1522 = sin(qJ(1,1));
t1528 = cos(qJ(1,1));
t1469 = t1507 * t1528 + t1510 * t1522;
t1521 = sin(qJ(2,1));
t1504 = t1521 * qJ(3,1);
t1513 = t1528 * pkin(4);
t1481 = t1522 * t1504 + t1513;
t1580 = t1522 * pkin(4);
t1482 = t1528 * t1504 - t1580;
t1527 = cos(qJ(2,1));
t1529 = pkin(1) + pkin(2);
t1545 = t1529 * t1527;
t1427 = t1469 * t1545 + t1481 * t1510 + t1507 * t1482;
t1468 = -t1507 * t1522 + t1510 * t1528;
t1488 = t1504 + t1545;
t1485 = 0.1e1 / t1488;
t1473 = t1507 * g(1) - t1510 * g(2);
t1476 = t1510 * g(1) + t1507 * g(2);
t1533 = t1473 * t1528 + t1476 * t1522;
t1570 = t1533 * t1521;
t1457 = (g(1) * t1522 - g(2) * t1528) * t1507 - (g(1) * t1528 + g(2) * t1522) * t1510;
t1439 = -g(3) * t1521 + t1457 * t1527;
t1532 = 0.1e1 / qJ(3,1);
t1548 = t1527 * t1532;
t1592 = t1439 * t1548;
t1605 = (t1427 * t1592 + t1468 * t1570) * t1485;
t1426 = t1468 * t1545 - t1507 * t1481 + t1482 * t1510;
t1604 = (t1426 * t1592 - t1469 * t1570) * t1485;
t1514 = legFrame(3,3);
t1505 = sin(t1514);
t1508 = cos(t1514);
t1518 = sin(qJ(1,3));
t1524 = cos(qJ(1,3));
t1464 = -t1505 * t1518 + t1508 * t1524;
t1517 = sin(qJ(2,3));
t1502 = t1517 * qJ(3,3);
t1511 = t1524 * pkin(4);
t1477 = t1518 * t1502 + t1511;
t1582 = t1518 * pkin(4);
t1478 = t1524 * t1502 - t1582;
t1523 = cos(qJ(2,3));
t1547 = t1529 * t1523;
t1422 = t1464 * t1547 - t1505 * t1477 + t1478 * t1508;
t1495 = t1518 * g(1) - t1524 * g(2);
t1496 = t1524 * g(1) + t1518 * g(2);
t1452 = t1495 * t1505 - t1496 * t1508;
t1431 = -g(3) * t1517 + t1452 * t1523;
t1530 = 0.1e1 / qJ(3,3);
t1550 = t1523 * t1530;
t1596 = t1431 * t1550;
t1603 = t1422 * t1596;
t1465 = t1505 * t1524 + t1508 * t1518;
t1423 = t1465 * t1547 + t1477 * t1508 + t1505 * t1478;
t1602 = t1423 * t1596;
t1515 = legFrame(2,3);
t1506 = sin(t1515);
t1509 = cos(t1515);
t1520 = sin(qJ(1,2));
t1526 = cos(qJ(1,2));
t1466 = -t1506 * t1520 + t1509 * t1526;
t1519 = sin(qJ(2,2));
t1503 = t1519 * qJ(3,2);
t1512 = t1526 * pkin(4);
t1479 = t1520 * t1503 + t1512;
t1581 = t1520 * pkin(4);
t1480 = t1526 * t1503 - t1581;
t1525 = cos(qJ(2,2));
t1546 = t1529 * t1525;
t1424 = t1466 * t1546 - t1506 * t1479 + t1480 * t1509;
t1497 = t1520 * g(1) - t1526 * g(2);
t1498 = t1526 * g(1) + t1520 * g(2);
t1455 = t1497 * t1506 - t1498 * t1509;
t1435 = -g(3) * t1519 + t1455 * t1525;
t1531 = 0.1e1 / qJ(3,2);
t1549 = t1525 * t1531;
t1594 = t1435 * t1549;
t1601 = t1424 * t1594;
t1467 = t1506 * t1526 + t1509 * t1520;
t1425 = t1467 * t1546 + t1479 * t1509 + t1506 * t1480;
t1600 = t1425 * t1594;
t1430 = g(3) * t1523 + t1452 * t1517;
t1597 = t1430 * t1530;
t1434 = g(3) * t1525 + t1455 * t1519;
t1595 = t1434 * t1531;
t1438 = g(3) * t1527 + t1457 * t1521;
t1593 = t1438 * t1532;
t1551 = t1521 * t1532;
t1552 = t1519 * t1531;
t1553 = t1517 * t1530;
t1591 = t1431 * t1553 + t1435 * t1552 + t1439 * t1551;
t1471 = t1505 * g(1) - t1508 * g(2);
t1474 = t1508 * g(1) + t1505 * g(2);
t1440 = t1518 * t1471 - t1474 * t1524;
t1472 = t1506 * g(1) - t1509 * g(2);
t1475 = t1509 * g(1) + t1506 * g(2);
t1441 = t1520 * t1472 - t1475 * t1526;
t1442 = t1522 * t1473 - t1476 * t1528;
t1557 = t1469 * t1485;
t1487 = t1503 + t1546;
t1484 = 0.1e1 / t1487;
t1560 = t1467 * t1484;
t1486 = t1502 + t1547;
t1483 = 0.1e1 / t1486;
t1564 = t1465 * t1483;
t1590 = t1440 * t1564 + t1441 * t1560 + t1442 * t1557;
t1558 = t1468 * t1485;
t1562 = t1466 * t1484;
t1566 = t1464 * t1483;
t1589 = t1440 * t1566 + t1441 * t1562 + t1442 * t1558;
t1579 = t1523 * qJ(3,3);
t1578 = t1525 * qJ(3,2);
t1577 = t1527 * qJ(3,1);
t1494 = t1527 * pkin(1) + t1504;
t1569 = t1533 * t1494;
t1450 = t1495 * t1508 + t1496 * t1505;
t1492 = t1523 * pkin(1) + t1502;
t1568 = t1450 * t1492;
t1453 = t1497 * t1509 + t1498 * t1506;
t1493 = t1525 * pkin(1) + t1503;
t1567 = t1453 * t1493;
t1565 = t1464 * t1517;
t1563 = t1465 * t1517;
t1561 = t1466 * t1519;
t1559 = t1467 * t1519;
t1556 = t1483 * t1523;
t1555 = t1484 * t1525;
t1554 = t1485 * t1527;
t1419 = -g(3) * t1492 - t1440 * (t1517 * pkin(1) - t1579);
t1544 = t1419 * t1550;
t1420 = -g(3) * t1493 - t1441 * (t1519 * pkin(1) - t1578);
t1543 = t1420 * t1549;
t1421 = -g(3) * t1494 - t1442 * (t1521 * pkin(1) - t1577);
t1542 = t1421 * t1548;
t1463 = t1488 * t1528 - t1580;
t1462 = t1487 * t1526 - t1581;
t1461 = t1486 * t1524 - t1582;
t1460 = t1488 * t1522 + t1513;
t1459 = t1487 * t1520 + t1512;
t1458 = t1486 * t1518 + t1511;
t1444 = t1472 * t1526 + t1475 * t1520;
t1443 = t1471 * t1524 + t1474 * t1518;
t1418 = -t1430 * t1553 - t1434 * t1552 - t1438 * t1551;
t1417 = (-t1426 * t1593 - t1469 * t1533) * t1554 + (-t1424 * t1595 - t1453 * t1467) * t1555 + (-t1422 * t1597 - t1450 * t1465) * t1556;
t1416 = (-t1427 * t1593 + t1468 * t1533) * t1554 + (-t1425 * t1595 + t1453 * t1466) * t1555 + (-t1423 * t1597 + t1450 * t1464) * t1556;
t1 = [0, -t1443 * t1564 - t1444 * t1560 - t1533 * t1557, t1590, 0, 0, 0, 0, 0, t1417, -t1604 + (t1444 * t1559 - t1601) * t1484 + (t1443 * t1563 - t1603) * t1483, t1417, -t1590, t1604 + (-t1453 * t1559 + t1601) * t1484 + (-t1450 * t1563 + t1603) * t1483, (-t1507 * t1460 + t1510 * t1463) * t1593 + (-t1506 * t1459 + t1509 * t1462) * t1595 + (-t1505 * t1458 + t1508 * t1461) * t1597 + (t1426 * t1542 - t1469 * t1569) * t1485 + (t1424 * t1543 - t1467 * t1567) * t1484 + (t1422 * t1544 - t1465 * t1568) * t1483, -g(1); 0, t1443 * t1566 + t1444 * t1562 + t1533 * t1558, -t1589, 0, 0, 0, 0, 0, t1416, -t1605 + (-t1444 * t1561 - t1600) * t1484 + (-t1443 * t1565 - t1602) * t1483, t1416, t1589, t1605 + (t1453 * t1561 + t1600) * t1484 + (t1450 * t1565 + t1602) * t1483, (t1460 * t1510 + t1507 * t1463) * t1593 + (t1459 * t1509 + t1506 * t1462) * t1595 + (t1458 * t1508 + t1505 * t1461) * t1597 + (t1427 * t1542 + t1468 * t1569) * t1485 + (t1425 * t1543 + t1466 * t1567) * t1484 + (t1423 * t1544 + t1464 * t1568) * t1483, -g(2); 0, 0, 0, 0, 0, 0, 0, 0, t1418, -t1591, t1418, 0, t1591, (t1521 * t1421 + (t1529 * t1521 - t1577) * t1438) * t1532 + (t1519 * t1420 + (t1529 * t1519 - t1578) * t1434) * t1531 + (t1517 * t1419 + (t1529 * t1517 - t1579) * t1430) * t1530, -g(3);];
tau_reg  = t1;
