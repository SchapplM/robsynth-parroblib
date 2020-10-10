% Calculate Gravitation load for parallel robot
% P4PRRRR8V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 22:59:48
% EndTime: 2020-09-20 22:59:49
% DurationCPUTime: 1.76s
% Computational Cost: add. (1315->201), mult. (2993->389), div. (120->9), fcn. (2654->30), ass. (0->189)
t1547 = sin(qJ(2,1));
t1553 = cos(qJ(2,1));
t1541 = legFrame(1,2);
t1517 = sin(t1541);
t1521 = cos(t1541);
t1498 = g(1) * t1521 - g(2) * t1517;
t1529 = sin(pkin(6));
t1531 = cos(pkin(6));
t1644 = g(3) * t1531;
t1480 = -t1498 * t1529 - t1644;
t1494 = g(1) * t1517 + g(2) * t1521;
t1530 = sin(pkin(3));
t1532 = cos(pkin(3));
t1586 = t1480 * t1532 + t1494 * t1530;
t1645 = g(3) * t1529;
t1590 = t1498 * t1531 - t1645;
t1566 = t1547 * t1586 + t1590 * t1553;
t1545 = sin(qJ(2,2));
t1551 = cos(qJ(2,2));
t1540 = legFrame(2,2);
t1516 = sin(t1540);
t1520 = cos(t1540);
t1497 = g(1) * t1520 - g(2) * t1516;
t1479 = -t1497 * t1529 - t1644;
t1493 = g(1) * t1516 + g(2) * t1520;
t1587 = t1479 * t1532 + t1493 * t1530;
t1591 = t1497 * t1531 - t1645;
t1567 = t1545 * t1587 + t1591 * t1551;
t1543 = sin(qJ(2,3));
t1549 = cos(qJ(2,3));
t1539 = legFrame(3,2);
t1515 = sin(t1539);
t1519 = cos(t1539);
t1496 = g(1) * t1519 - g(2) * t1515;
t1478 = -t1496 * t1529 - t1644;
t1492 = g(1) * t1515 + g(2) * t1519;
t1588 = t1478 * t1532 + t1492 * t1530;
t1592 = t1496 * t1531 - t1645;
t1568 = t1543 * t1588 + t1592 * t1549;
t1534 = sin(qJ(2,4));
t1536 = cos(qJ(2,4));
t1538 = legFrame(4,2);
t1514 = sin(t1538);
t1518 = cos(t1538);
t1495 = g(1) * t1518 - g(2) * t1514;
t1477 = -t1495 * t1529 - t1644;
t1491 = g(1) * t1514 + g(2) * t1518;
t1589 = t1477 * t1532 + t1491 * t1530;
t1593 = t1495 * t1531 - t1645;
t1569 = t1534 * t1589 + t1593 * t1536;
t1535 = cos(qJ(3,4));
t1652 = pkin(2) * t1535;
t1489 = -t1536 * pkin(5) + t1534 * t1652;
t1533 = sin(qJ(3,4));
t1653 = pkin(2) * t1533;
t1470 = t1489 * t1530 + t1532 * t1653;
t1657 = 0.1e1 / t1470;
t1548 = cos(qJ(3,3));
t1648 = pkin(2) * t1548;
t1499 = -t1549 * pkin(5) + t1543 * t1648;
t1542 = sin(qJ(3,3));
t1651 = pkin(2) * t1542;
t1474 = t1499 * t1530 + t1532 * t1651;
t1656 = 0.1e1 / t1474;
t1550 = cos(qJ(3,2));
t1647 = pkin(2) * t1550;
t1500 = -t1551 * pkin(5) + t1545 * t1647;
t1544 = sin(qJ(3,2));
t1650 = pkin(2) * t1544;
t1475 = t1500 * t1530 + t1532 * t1650;
t1655 = 0.1e1 / t1475;
t1552 = cos(qJ(3,1));
t1646 = pkin(2) * t1552;
t1501 = -t1553 * pkin(5) + t1547 * t1646;
t1546 = sin(qJ(3,1));
t1649 = pkin(2) * t1546;
t1476 = t1501 * t1530 + t1532 * t1649;
t1654 = 0.1e1 / t1476;
t1643 = t1657 / t1535;
t1642 = t1656 / t1548;
t1641 = t1655 / t1550;
t1640 = t1654 / t1552;
t1639 = t1657 * t1491;
t1638 = t1656 * t1492;
t1637 = t1655 * t1493;
t1636 = t1654 * t1494;
t1630 = t1491 * t1532;
t1628 = t1492 * t1532;
t1626 = t1493 * t1532;
t1624 = t1494 * t1532;
t1623 = mrSges(3,2) * t1530 * t1644;
t1622 = t1529 * t1530;
t1621 = t1532 * t1534;
t1620 = t1532 * t1536;
t1619 = t1532 * t1543;
t1618 = t1532 * t1545;
t1617 = t1532 * t1547;
t1616 = t1532 * t1549;
t1615 = t1532 * t1551;
t1614 = t1532 * t1553;
t1613 = t1533 * t1536;
t1612 = t1542 * t1549;
t1611 = t1544 * t1551;
t1610 = t1546 * t1553;
t1437 = ((t1477 * t1530 - t1630) * mrSges(3,1) + t1569 * mrSges(3,2)) * t1535 + t1533 * (t1623 + (t1495 * t1622 + t1630) * mrSges(3,2) + t1569 * mrSges(3,1));
t1609 = t1437 * t1643;
t1438 = ((t1478 * t1530 - t1628) * mrSges(3,1) + t1568 * mrSges(3,2)) * t1548 + t1542 * (t1623 + (t1496 * t1622 + t1628) * mrSges(3,2) + t1568 * mrSges(3,1));
t1608 = t1438 * t1642;
t1439 = ((t1479 * t1530 - t1626) * mrSges(3,1) + t1567 * mrSges(3,2)) * t1550 + t1544 * (t1623 + (t1497 * t1622 + t1626) * mrSges(3,2) + t1567 * mrSges(3,1));
t1607 = t1439 * t1641;
t1440 = ((t1480 * t1530 - t1624) * mrSges(3,1) + t1566 * mrSges(3,2)) * t1552 + t1546 * (t1623 + (t1498 * t1622 + t1624) * mrSges(3,2) + t1566 * mrSges(3,1));
t1606 = t1440 * t1640;
t1537 = mrSges(2,2) - mrSges(3,3);
t1441 = t1569 * t1537 + (t1534 * t1593 - t1589 * t1536) * (mrSges(3,1) * t1535 - mrSges(3,2) * t1533 + mrSges(2,1));
t1605 = t1441 * t1643;
t1442 = t1568 * t1537 + (t1543 * t1592 - t1588 * t1549) * (mrSges(3,1) * t1548 - mrSges(3,2) * t1542 + mrSges(2,1));
t1604 = t1442 * t1642;
t1443 = t1567 * t1537 + (t1545 * t1591 - t1587 * t1551) * (mrSges(3,1) * t1550 - mrSges(3,2) * t1544 + mrSges(2,1));
t1603 = t1443 * t1641;
t1444 = t1566 * t1537 + (t1547 * t1590 - t1586 * t1553) * (mrSges(3,1) * t1552 - mrSges(3,2) * t1546 + mrSges(2,1));
t1602 = t1444 * t1640;
t1457 = (-t1529 * t1534 + t1531 * t1620) * t1652 + pkin(5) * (t1529 * t1536 + t1531 * t1621);
t1601 = t1457 * t1609;
t1458 = (-t1529 * t1543 + t1531 * t1616) * t1648 + pkin(5) * (t1529 * t1549 + t1531 * t1619);
t1600 = t1458 * t1608;
t1459 = (-t1529 * t1545 + t1531 * t1615) * t1647 + pkin(5) * (t1529 * t1551 + t1531 * t1618);
t1599 = t1459 * t1607;
t1460 = (-t1529 * t1547 + t1531 * t1614) * t1646 + pkin(5) * (t1529 * t1553 + t1531 * t1617);
t1598 = t1460 * t1606;
t1581 = t1530 * t1535 + t1533 * t1621;
t1461 = t1529 * t1613 + t1581 * t1531;
t1597 = t1461 * t1605;
t1580 = t1530 * t1548 + t1542 * t1619;
t1462 = t1529 * t1612 + t1580 * t1531;
t1596 = t1462 * t1604;
t1579 = t1530 * t1550 + t1544 * t1618;
t1463 = t1529 * t1611 + t1579 * t1531;
t1595 = t1463 * t1603;
t1578 = t1530 * t1552 + t1546 * t1617;
t1464 = t1529 * t1610 + t1578 * t1531;
t1594 = t1464 * t1602;
t1585 = -t1489 * t1532 + t1530 * t1653;
t1584 = -t1499 * t1532 + t1530 * t1651;
t1583 = -t1500 * t1532 + t1530 * t1650;
t1582 = -t1501 * t1532 + t1530 * t1649;
t1554 = xP(4);
t1522 = sin(t1554);
t1523 = cos(t1554);
t1557 = koppelP(4,2);
t1561 = koppelP(4,1);
t1481 = -t1522 * t1561 - t1523 * t1557;
t1485 = -t1522 * t1557 + t1523 * t1561;
t1573 = (-t1481 * t1518 + t1485 * t1514) * t1643;
t1558 = koppelP(3,2);
t1562 = koppelP(3,1);
t1482 = -t1522 * t1562 - t1523 * t1558;
t1486 = -t1522 * t1558 + t1523 * t1562;
t1572 = (-t1482 * t1519 + t1486 * t1515) * t1642;
t1559 = koppelP(2,2);
t1563 = koppelP(2,1);
t1483 = -t1522 * t1563 - t1523 * t1559;
t1487 = -t1522 * t1559 + t1523 * t1563;
t1571 = (-t1483 * t1520 + t1487 * t1516) * t1641;
t1560 = koppelP(1,2);
t1564 = koppelP(1,1);
t1484 = -t1522 * t1564 - t1523 * t1560;
t1488 = -t1522 * t1560 + t1523 * t1564;
t1570 = (-t1484 * t1521 + t1488 * t1517) * t1640;
t1565 = 0.1e1 / pkin(2);
t1556 = mrSges(4,1);
t1555 = mrSges(4,2);
t1525 = m(1) + m(2) + m(3);
t1504 = pkin(5) * t1547 + t1553 * t1646;
t1503 = pkin(5) * t1545 + t1551 * t1647;
t1502 = pkin(5) * t1543 + t1549 * t1648;
t1490 = pkin(5) * t1534 + t1536 * t1652;
t1456 = t1504 * t1531 + t1582 * t1529;
t1455 = t1503 * t1531 + t1583 * t1529;
t1454 = t1502 * t1531 + t1584 * t1529;
t1453 = t1490 * t1531 + t1585 * t1529;
t1452 = -t1456 * t1517 + t1476 * t1521;
t1451 = t1456 * t1521 + t1476 * t1517;
t1450 = -t1455 * t1516 + t1475 * t1520;
t1449 = t1455 * t1520 + t1475 * t1516;
t1448 = -t1454 * t1515 + t1474 * t1519;
t1447 = t1454 * t1519 + t1474 * t1515;
t1446 = -t1453 * t1514 + t1470 * t1518;
t1445 = t1453 * t1518 + t1470 * t1514;
t1 = [-t1518 * t1597 - t1519 * t1596 - t1520 * t1595 - t1521 * t1594 - g(1) * m(4) + (-t1518 * t1601 - t1519 * t1600 - t1520 * t1599 - t1521 * t1598) * t1565 + (-t1445 * t1639 - t1447 * t1638 - t1449 * t1637 - t1451 * t1636) * t1525; t1514 * t1597 + t1515 * t1596 + t1516 * t1595 + t1517 * t1594 - g(2) * m(4) + (t1514 * t1601 + t1515 * t1600 + t1516 * t1599 + t1517 * t1598) * t1565 + (-t1446 * t1639 - t1448 * t1638 - t1450 * t1637 - t1452 * t1636) * t1525; (t1578 * t1529 - t1531 * t1610) * t1602 + (t1579 * t1529 - t1531 * t1611) * t1603 + (t1580 * t1529 - t1531 * t1612) * t1604 + (t1581 * t1529 - t1531 * t1613) * t1605 - g(3) * m(4) + (((t1529 * t1614 + t1531 * t1547) * t1646 + (t1529 * t1617 - t1531 * t1553) * pkin(5)) * t1606 + ((t1529 * t1615 + t1531 * t1545) * t1647 + (t1529 * t1618 - t1531 * t1551) * pkin(5)) * t1607 + ((t1529 * t1616 + t1531 * t1543) * t1648 + (t1529 * t1619 - t1531 * t1549) * pkin(5)) * t1608 + ((t1529 * t1620 + t1531 * t1534) * t1652 + (t1529 * t1621 - t1531 * t1536) * pkin(5)) * t1609) * t1565 + (-(-t1504 * t1529 + t1582 * t1531) * t1636 - (-t1503 * t1529 + t1583 * t1531) * t1637 - (-t1529 * t1502 + t1584 * t1531) * t1638 - (-t1490 * t1529 + t1585 * t1531) * t1639) * t1525; -(-g(1) * t1556 - g(2) * t1555) * t1522 + t1523 * (g(1) * t1555 - g(2) * t1556) + t1444 * t1464 * t1570 + t1443 * t1463 * t1571 + t1442 * t1462 * t1572 + t1441 * t1461 * t1573 + (-(t1451 * t1484 + t1452 * t1488) * t1636 - (t1449 * t1483 + t1450 * t1487) * t1637 - (t1447 * t1482 + t1448 * t1486) * t1638 - (t1445 * t1481 + t1446 * t1485) * t1639) * t1525 + (t1437 * t1457 * t1573 + t1438 * t1458 * t1572 + t1439 * t1459 * t1571 + t1460 * t1440 * t1570) * t1565;];
taugX  = t1;
