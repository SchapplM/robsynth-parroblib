% Calculate Gravitation load for parallel robot
% P4PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:10:21
% EndTime: 2020-08-07 11:10:23
% DurationCPUTime: 2.26s
% Computational Cost: add. (1638->245), mult. (3131->452), div. (64->9), fcn. (2810->30), ass. (0->209)
t1573 = sin(qJ(2,4));
t1575 = cos(qJ(2,4));
t1568 = legFrame(4,3);
t1549 = sin(t1568);
t1553 = cos(t1568);
t1524 = -t1549 * g(1) + t1553 * g(2);
t1528 = t1553 * g(1) + t1549 * g(2);
t1564 = sin(pkin(8));
t1566 = cos(pkin(8));
t1617 = t1524 * t1564 + t1528 * t1566;
t1567 = cos(pkin(4));
t1618 = t1524 * t1566 - t1528 * t1564;
t1565 = sin(pkin(4));
t1667 = g(3) * t1565;
t1674 = t1618 * t1567 + t1667;
t1606 = t1674 * t1573 + t1617 * t1575;
t1577 = sin(qJ(2,3));
t1583 = cos(qJ(2,3));
t1569 = legFrame(3,3);
t1550 = sin(t1569);
t1554 = cos(t1569);
t1525 = -t1550 * g(1) + t1554 * g(2);
t1529 = t1554 * g(1) + t1550 * g(2);
t1615 = t1525 * t1564 + t1529 * t1566;
t1616 = t1525 * t1566 - t1529 * t1564;
t1675 = t1616 * t1567 + t1667;
t1605 = t1675 * t1577 + t1615 * t1583;
t1579 = sin(qJ(2,2));
t1585 = cos(qJ(2,2));
t1570 = legFrame(2,3);
t1551 = sin(t1570);
t1555 = cos(t1570);
t1526 = -t1551 * g(1) + t1555 * g(2);
t1530 = t1555 * g(1) + t1551 * g(2);
t1613 = t1526 * t1564 + t1530 * t1566;
t1614 = t1526 * t1566 - t1530 * t1564;
t1676 = t1614 * t1567 + t1667;
t1604 = t1676 * t1579 + t1613 * t1585;
t1581 = sin(qJ(2,1));
t1587 = cos(qJ(2,1));
t1571 = legFrame(1,3);
t1552 = sin(t1571);
t1556 = cos(t1571);
t1527 = -t1552 * g(1) + t1556 * g(2);
t1531 = t1556 * g(1) + t1552 * g(2);
t1611 = t1527 * t1564 + t1531 * t1566;
t1612 = t1527 * t1566 - t1531 * t1564;
t1677 = t1612 * t1567 + t1667;
t1603 = t1677 * t1581 + t1611 * t1587;
t1666 = g(3) * t1567;
t1589 = pkin(7) + pkin(6);
t1628 = t1581 * t1589;
t1539 = pkin(2) * t1587 + t1628;
t1547 = t1589 * t1587;
t1536 = pkin(2) * t1581 - t1547;
t1580 = sin(qJ(3,1));
t1643 = t1565 * t1580;
t1607 = pkin(3) * t1643 - t1536 * t1567;
t1681 = t1539 * t1566 + t1607 * t1564;
t1629 = t1579 * t1589;
t1538 = pkin(2) * t1585 + t1629;
t1546 = t1589 * t1585;
t1535 = pkin(2) * t1579 - t1546;
t1578 = sin(qJ(3,2));
t1644 = t1565 * t1578;
t1608 = pkin(3) * t1644 - t1535 * t1567;
t1680 = t1538 * t1566 + t1608 * t1564;
t1630 = t1577 * t1589;
t1537 = pkin(2) * t1583 + t1630;
t1545 = t1589 * t1583;
t1534 = pkin(2) * t1577 - t1545;
t1576 = sin(qJ(3,3));
t1645 = t1565 * t1576;
t1609 = pkin(3) * t1645 - t1534 * t1567;
t1679 = t1537 * t1566 + t1609 * t1564;
t1631 = t1573 * t1589;
t1533 = pkin(2) * t1575 + t1631;
t1541 = t1589 * t1575;
t1532 = pkin(2) * t1573 - t1541;
t1572 = sin(qJ(3,4));
t1647 = t1565 * t1572;
t1610 = pkin(3) * t1647 - t1532 * t1567;
t1678 = t1533 * t1566 + t1610 * t1564;
t1574 = cos(qJ(3,4));
t1559 = t1574 ^ 2;
t1673 = pkin(3) * t1559;
t1582 = cos(qJ(3,3));
t1561 = t1582 ^ 2;
t1672 = pkin(3) * t1561;
t1584 = cos(qJ(3,2));
t1562 = t1584 ^ 2;
t1671 = pkin(3) * t1562;
t1586 = cos(qJ(3,1));
t1563 = t1586 ^ 2;
t1670 = pkin(3) * t1563;
t1669 = pkin(3) * t1565;
t1560 = m(1) + m(2) + m(3);
t1668 = g(3) * t1560;
t1665 = mrSges(3,2) * t1565;
t1664 = m(3) * pkin(2) + mrSges(2,1);
t1540 = t1574 * pkin(3) + pkin(2);
t1505 = t1573 * t1540 - t1541;
t1627 = mrSges(3,2) * t1666;
t1639 = t1567 * t1572;
t1646 = t1565 * t1574;
t1663 = (((t1618 * t1565 - t1666) * mrSges(3,1) + t1606 * mrSges(3,2)) * t1574 + (t1606 * mrSges(3,1) - t1618 * t1665 + t1627) * t1572) / (t1505 * t1646 + t1540 * t1639);
t1542 = t1582 * pkin(3) + pkin(2);
t1509 = t1577 * t1542 - t1545;
t1637 = t1567 * t1576;
t1642 = t1565 * t1582;
t1662 = (((t1616 * t1565 - t1666) * mrSges(3,1) + t1605 * mrSges(3,2)) * t1582 + (t1605 * mrSges(3,1) - t1616 * t1665 + t1627) * t1576) / (t1509 * t1642 + t1542 * t1637);
t1543 = t1584 * pkin(3) + pkin(2);
t1510 = t1579 * t1543 - t1546;
t1635 = t1567 * t1578;
t1641 = t1565 * t1584;
t1661 = (((t1614 * t1565 - t1666) * mrSges(3,1) + t1604 * mrSges(3,2)) * t1584 + (t1604 * mrSges(3,1) - t1614 * t1665 + t1627) * t1578) / (t1510 * t1641 + t1543 * t1635);
t1544 = t1586 * pkin(3) + pkin(2);
t1511 = t1581 * t1544 - t1547;
t1633 = t1567 * t1580;
t1640 = t1565 * t1586;
t1660 = (((t1612 * t1565 - t1666) * mrSges(3,1) + t1603 * mrSges(3,2)) * t1586 + (t1603 * mrSges(3,1) - t1612 * t1665 + t1627) * t1580) / (t1511 * t1640 + t1544 * t1633);
t1472 = 0.1e1 / (t1573 * t1559 * t1669 + (pkin(3) * t1639 + t1532 * t1565) * t1574 + pkin(2) * t1639);
t1548 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1659 = (-t1606 * t1548 + (t1617 * t1573 - t1575 * t1674) * (t1574 * mrSges(3,1) - mrSges(3,2) * t1572 + t1664)) * t1472;
t1473 = 0.1e1 / (t1577 * t1561 * t1669 + (pkin(3) * t1637 + t1534 * t1565) * t1582 + pkin(2) * t1637);
t1658 = (-t1605 * t1548 + (t1615 * t1577 - t1583 * t1675) * (t1582 * mrSges(3,1) - mrSges(3,2) * t1576 + t1664)) * t1473;
t1474 = 0.1e1 / (t1579 * t1562 * t1669 + (pkin(3) * t1635 + t1535 * t1565) * t1584 + pkin(2) * t1635);
t1657 = (-t1604 * t1548 + (t1613 * t1579 - t1585 * t1676) * (t1584 * mrSges(3,1) - mrSges(3,2) * t1578 + t1664)) * t1474;
t1475 = 0.1e1 / (t1581 * t1563 * t1669 + (pkin(3) * t1633 + t1536 * t1565) * t1586 + pkin(2) * t1633);
t1656 = (-t1603 * t1548 + (t1611 * t1581 - t1587 * t1677) * (t1586 * mrSges(3,1) - mrSges(3,2) * t1580 + t1664)) * t1475;
t1655 = (t1540 * t1575 + t1631) * t1567;
t1654 = (t1542 * t1583 + t1630) * t1567;
t1653 = (t1543 * t1585 + t1629) * t1567;
t1652 = (t1544 * t1587 + t1628) * t1567;
t1638 = t1567 * t1573;
t1636 = t1567 * t1577;
t1634 = t1567 * t1579;
t1632 = t1567 * t1581;
t1626 = pkin(2) * t1647;
t1625 = pkin(2) * t1645;
t1624 = pkin(2) * t1644;
t1623 = pkin(2) * t1643;
t1602 = 0.1e1 / pkin(3);
t1601 = koppelP(1,1);
t1600 = koppelP(2,1);
t1599 = koppelP(3,1);
t1598 = koppelP(4,1);
t1597 = koppelP(1,2);
t1596 = koppelP(2,2);
t1595 = koppelP(3,2);
t1594 = koppelP(4,2);
t1593 = mrSges(4,1);
t1592 = mrSges(4,2);
t1591 = xP(4);
t1558 = cos(t1591);
t1557 = sin(t1591);
t1519 = -t1557 * t1597 + t1558 * t1601;
t1518 = -t1557 * t1596 + t1558 * t1600;
t1517 = -t1557 * t1595 + t1558 * t1599;
t1516 = -t1557 * t1594 + t1558 * t1598;
t1515 = -t1557 * t1601 - t1558 * t1597;
t1514 = -t1557 * t1600 - t1558 * t1596;
t1513 = -t1557 * t1599 - t1558 * t1595;
t1512 = -t1557 * t1598 - t1558 * t1594;
t1503 = t1564 * t1587 + t1566 * t1632;
t1502 = t1564 * t1585 + t1566 * t1634;
t1501 = t1564 * t1583 + t1566 * t1636;
t1500 = t1564 * t1632 - t1566 * t1587;
t1499 = t1564 * t1634 - t1566 * t1585;
t1498 = t1564 * t1636 - t1566 * t1583;
t1497 = t1564 * t1575 + t1566 * t1638;
t1496 = t1564 * t1638 - t1566 * t1575;
t1495 = t1566 * t1552 + t1556 * t1564;
t1494 = t1566 * t1551 + t1555 * t1564;
t1493 = t1566 * t1550 + t1554 * t1564;
t1492 = t1566 * t1549 + t1553 * t1564;
t1491 = -t1564 * t1552 + t1556 * t1566;
t1490 = -t1564 * t1551 + t1555 * t1566;
t1489 = -t1564 * t1550 + t1554 * t1566;
t1488 = -t1564 * t1549 + t1553 * t1566;
t1479 = t1564 * t1539 - t1607 * t1566;
t1478 = t1564 * t1538 - t1608 * t1566;
t1477 = t1564 * t1537 - t1609 * t1566;
t1476 = t1564 * t1533 - t1610 * t1566;
t1471 = -t1511 * t1491 - t1495 * t1652;
t1470 = -t1510 * t1490 - t1494 * t1653;
t1469 = -t1509 * t1489 - t1493 * t1654;
t1468 = -t1491 * t1652 + t1511 * t1495;
t1467 = -t1490 * t1653 + t1510 * t1494;
t1466 = -t1489 * t1654 + t1509 * t1493;
t1465 = -t1505 * t1488 - t1492 * t1655;
t1464 = -t1488 * t1655 + t1505 * t1492;
t1463 = -t1495 * t1640 - (-t1587 * t1491 + t1495 * t1632) * t1580;
t1462 = -t1494 * t1641 - (-t1585 * t1490 + t1494 * t1634) * t1578;
t1461 = -t1493 * t1642 - (-t1583 * t1489 + t1493 * t1636) * t1576;
t1460 = -t1491 * t1640 - (t1491 * t1632 + t1587 * t1495) * t1580;
t1459 = -t1490 * t1641 - (t1490 * t1634 + t1585 * t1494) * t1578;
t1458 = -t1489 * t1642 - (t1489 * t1636 + t1583 * t1493) * t1576;
t1457 = -t1492 * t1646 - (-t1575 * t1488 + t1492 * t1638) * t1572;
t1456 = -t1488 * t1646 - (t1488 * t1638 + t1575 * t1492) * t1572;
t1455 = (-t1552 * t1500 + t1503 * t1556) * t1670 + (t1479 * t1556 + t1681 * t1552) * t1586 - t1491 * t1623;
t1454 = (-t1551 * t1499 + t1502 * t1555) * t1671 + (t1478 * t1555 + t1680 * t1551) * t1584 - t1490 * t1624;
t1453 = (-t1550 * t1498 + t1501 * t1554) * t1672 + (t1477 * t1554 + t1679 * t1550) * t1582 - t1489 * t1625;
t1452 = -(t1500 * t1556 + t1552 * t1503) * t1670 + (-t1552 * t1479 + t1681 * t1556) * t1586 + t1495 * t1623;
t1451 = -(t1499 * t1555 + t1551 * t1502) * t1671 + (-t1551 * t1478 + t1680 * t1555) * t1584 + t1494 * t1624;
t1450 = -(t1498 * t1554 + t1550 * t1501) * t1672 + (-t1550 * t1477 + t1679 * t1554) * t1582 + t1493 * t1625;
t1449 = (-t1549 * t1496 + t1497 * t1553) * t1673 + (t1476 * t1553 + t1678 * t1549) * t1574 - t1488 * t1626;
t1448 = -(t1496 * t1553 + t1549 * t1497) * t1673 + (-t1549 * t1476 + t1678 * t1553) * t1574 + t1492 * t1626;
t1 = [t1456 * t1659 + t1458 * t1658 + t1459 * t1657 + t1460 * t1656 - g(1) * m(4) + (t1464 * t1663 + t1466 * t1662 + t1467 * t1661 + t1468 * t1660) * t1602 + (-t1448 * t1472 - t1450 * t1473 - t1451 * t1474 - t1452 * t1475) * t1668; t1457 * t1659 + t1461 * t1658 + t1462 * t1657 + t1463 * t1656 - g(2) * m(4) + (t1465 * t1663 + t1469 * t1662 + t1470 * t1661 + t1471 * t1660) * t1602 + (-t1449 * t1472 - t1453 * t1473 - t1454 * t1474 - t1455 * t1475) * t1668; (-m(4) - 0.4e1 * t1560) * g(3); -(-g(1) * t1593 - g(2) * t1592) * t1557 + t1558 * (g(1) * t1592 - g(2) * t1593) + (t1460 * t1515 + t1463 * t1519) * t1656 + (t1459 * t1514 + t1462 * t1518) * t1657 + (t1458 * t1513 + t1461 * t1517) * t1658 + (t1456 * t1512 + t1457 * t1516) * t1659 + ((-t1452 * t1515 - t1455 * t1519) * t1475 + (-t1451 * t1514 - t1454 * t1518) * t1474 + (-t1450 * t1513 - t1453 * t1517) * t1473 + (-t1448 * t1512 - t1449 * t1516) * t1472) * t1668 + ((t1468 * t1515 + t1471 * t1519) * t1660 + (t1467 * t1514 + t1470 * t1518) * t1661 + (t1466 * t1513 + t1469 * t1517) * t1662 + (t1464 * t1512 + t1465 * t1516) * t1663) * t1602;];
taugX  = t1;
