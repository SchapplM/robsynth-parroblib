% Calculate Gravitation load for parallel robot
% P4PRRRR8V2G3A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:25:28
% EndTime: 2020-08-07 11:25:31
% DurationCPUTime: 2.92s
% Computational Cost: add. (1851->256), mult. (3978->498), div. (80->5), fcn. (3150->30), ass. (0->201)
t1588 = sin(qJ(2,4));
t1590 = cos(qJ(2,4));
t1591 = legFrame(4,2);
t1568 = sin(t1591);
t1572 = cos(t1591);
t1553 = t1572 * g(1) - t1568 * g(2);
t1585 = cos(pkin(8));
t1583 = sin(pkin(8));
t1715 = g(3) * t1583;
t1648 = -t1553 * t1585 + t1715;
t1549 = t1568 * g(1) + t1572 * g(2);
t1584 = sin(pkin(4));
t1586 = cos(pkin(4));
t1567 = g(3) * t1585;
t1649 = -t1553 * t1583 - t1567;
t1733 = t1549 * t1584 + t1649 * t1586;
t1625 = t1733 * t1588 - t1648 * t1590;
t1596 = sin(qJ(2,3));
t1602 = cos(qJ(2,3));
t1592 = legFrame(3,2);
t1569 = sin(t1592);
t1573 = cos(t1592);
t1554 = t1573 * g(1) - t1569 * g(2);
t1646 = -t1554 * t1585 + t1715;
t1550 = t1569 * g(1) + t1573 * g(2);
t1647 = -t1554 * t1583 - t1567;
t1734 = t1550 * t1584 + t1647 * t1586;
t1624 = t1734 * t1596 - t1646 * t1602;
t1598 = sin(qJ(2,2));
t1604 = cos(qJ(2,2));
t1593 = legFrame(2,2);
t1570 = sin(t1593);
t1574 = cos(t1593);
t1555 = t1574 * g(1) - t1570 * g(2);
t1644 = -t1555 * t1585 + t1715;
t1551 = t1570 * g(1) + t1574 * g(2);
t1645 = -t1555 * t1583 - t1567;
t1735 = t1551 * t1584 + t1645 * t1586;
t1623 = t1735 * t1598 - t1644 * t1604;
t1600 = sin(qJ(2,1));
t1606 = cos(qJ(2,1));
t1594 = legFrame(1,2);
t1571 = sin(t1594);
t1575 = cos(t1594);
t1556 = t1575 * g(1) - t1571 * g(2);
t1642 = -t1556 * t1585 + t1715;
t1552 = t1571 * g(1) + t1575 * g(2);
t1643 = -t1556 * t1583 - t1567;
t1736 = t1552 * t1584 + t1643 * t1586;
t1622 = t1736 * t1600 - t1642 * t1606;
t1728 = m(3) / pkin(3);
t1587 = sin(qJ(3,4));
t1727 = pkin(2) * t1587;
t1595 = sin(qJ(3,3));
t1726 = pkin(2) * t1595;
t1597 = sin(qJ(3,2));
t1725 = pkin(2) * t1597;
t1599 = sin(qJ(3,1));
t1724 = pkin(2) * t1599;
t1589 = cos(qJ(3,4));
t1723 = pkin(3) * t1589 ^ 2;
t1601 = cos(qJ(3,3));
t1722 = pkin(3) * t1601 ^ 2;
t1603 = cos(qJ(3,2));
t1721 = pkin(3) * t1603 ^ 2;
t1605 = cos(qJ(3,1));
t1720 = pkin(3) * t1605 ^ 2;
t1719 = pkin(3) * t1589;
t1718 = pkin(3) * t1601;
t1717 = pkin(3) * t1603;
t1716 = pkin(3) * t1605;
t1608 = pkin(7) + pkin(6);
t1557 = pkin(2) * t1588 - t1608 * t1590;
t1670 = t1586 * t1587;
t1521 = pkin(3) * t1670 + t1557 * t1584;
t1682 = t1584 * t1588;
t1509 = 0.1e1 / (pkin(2) * t1670 + t1521 * t1589 + t1682 * t1723);
t1685 = t1583 * t1584;
t1686 = rSges(3,2) * t1584 * t1567;
t1700 = t1549 * t1586;
t1714 = (((t1649 * t1584 - t1700) * rSges(3,1) + t1625 * rSges(3,2)) * t1589 + (t1686 + (t1553 * t1685 + t1700) * rSges(3,2) + t1625 * rSges(3,1)) * t1587) * t1509;
t1559 = pkin(2) * t1596 - t1608 * t1602;
t1667 = t1586 * t1595;
t1522 = pkin(3) * t1667 + t1559 * t1584;
t1679 = t1584 * t1596;
t1510 = 0.1e1 / (pkin(2) * t1667 + t1522 * t1601 + t1679 * t1722);
t1697 = t1550 * t1586;
t1713 = (((t1647 * t1584 - t1697) * rSges(3,1) + t1624 * rSges(3,2)) * t1601 + (t1686 + (t1554 * t1685 + t1697) * rSges(3,2) + t1624 * rSges(3,1)) * t1595) * t1510;
t1560 = pkin(2) * t1598 - t1608 * t1604;
t1665 = t1586 * t1597;
t1523 = pkin(3) * t1665 + t1560 * t1584;
t1677 = t1584 * t1598;
t1511 = 0.1e1 / (pkin(2) * t1665 + t1523 * t1603 + t1677 * t1721);
t1694 = t1551 * t1586;
t1712 = (((t1645 * t1584 - t1694) * rSges(3,1) + t1623 * rSges(3,2)) * t1603 + (t1686 + (t1555 * t1685 + t1694) * rSges(3,2) + t1623 * rSges(3,1)) * t1597) * t1511;
t1561 = pkin(2) * t1600 - t1608 * t1606;
t1663 = t1586 * t1599;
t1524 = pkin(3) * t1663 + t1561 * t1584;
t1675 = t1584 * t1600;
t1512 = 0.1e1 / (pkin(2) * t1663 + t1524 * t1605 + t1675 * t1720);
t1691 = t1552 * t1586;
t1711 = (((t1643 * t1584 - t1691) * rSges(3,1) + t1622 * rSges(3,2)) * t1605 + (t1686 + (t1556 * t1685 + t1691) * rSges(3,2) + t1622 * rSges(3,1)) * t1599) * t1512;
t1565 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t1658 = m(2) * rSges(2,1) + pkin(2) * m(3);
t1493 = t1625 * t1565 + (-t1588 * t1648 - t1590 * t1733) * ((rSges(3,1) * t1589 - rSges(3,2) * t1587) * m(3) + t1658);
t1669 = t1586 * t1588;
t1531 = t1583 * t1590 + t1585 * t1669;
t1681 = t1584 * t1589;
t1710 = t1493 * (t1587 * t1531 + t1585 * t1681);
t1494 = t1624 * t1565 + (-t1596 * t1646 - t1602 * t1734) * ((rSges(3,1) * t1601 - rSges(3,2) * t1595) * m(3) + t1658);
t1666 = t1586 * t1596;
t1538 = t1583 * t1602 + t1585 * t1666;
t1674 = t1584 * t1601;
t1709 = t1494 * (t1595 * t1538 + t1585 * t1674);
t1495 = t1623 * t1565 + (-t1598 * t1644 - t1604 * t1735) * ((rSges(3,1) * t1603 - rSges(3,2) * t1597) * m(3) + t1658);
t1664 = t1586 * t1598;
t1539 = t1583 * t1604 + t1585 * t1664;
t1673 = t1584 * t1603;
t1708 = t1495 * (t1597 * t1539 + t1585 * t1673);
t1496 = t1622 * t1565 + (-t1600 * t1642 - t1606 * t1736) * ((rSges(3,1) * t1605 - rSges(3,2) * t1599) * m(3) + t1658);
t1662 = t1586 * t1600;
t1540 = t1583 * t1606 + t1585 * t1662;
t1672 = t1584 * t1605;
t1707 = t1496 * (t1599 * t1540 + t1585 * t1672);
t1706 = t1509 * t1549;
t1705 = t1510 * t1550;
t1704 = t1511 * t1551;
t1703 = t1512 * t1552;
t1579 = m(1) + m(2) + m(3);
t1702 = t1549 * t1579;
t1699 = t1550 * t1579;
t1696 = t1551 * t1579;
t1693 = t1552 * t1579;
t1684 = t1583 * t1586;
t1683 = t1584 * t1587;
t1680 = t1584 * t1595;
t1678 = t1584 * t1597;
t1676 = t1584 * t1599;
t1671 = t1585 * t1586;
t1668 = t1586 * t1590;
t1661 = t1586 * t1602;
t1660 = t1586 * t1604;
t1659 = t1586 * t1606;
t1558 = pkin(2) * t1590 + t1588 * t1608;
t1657 = ((t1583 * t1588 - t1585 * t1668) * t1719 - t1558 * t1671 + t1557 * t1583) * t1714;
t1562 = pkin(2) * t1602 + t1596 * t1608;
t1656 = ((t1583 * t1596 - t1585 * t1661) * t1718 - t1562 * t1671 + t1559 * t1583) * t1713;
t1563 = pkin(2) * t1604 + t1598 * t1608;
t1655 = ((t1583 * t1598 - t1585 * t1660) * t1717 - t1563 * t1671 + t1560 * t1583) * t1712;
t1564 = pkin(2) * t1606 + t1600 * t1608;
t1654 = ((t1583 * t1600 - t1585 * t1659) * t1716 - t1564 * t1671 + t1561 * t1583) * t1711;
t1653 = t1509 * t1710;
t1652 = t1510 * t1709;
t1651 = t1511 * t1708;
t1650 = t1512 * t1707;
t1610 = xP(4);
t1576 = sin(t1610);
t1577 = cos(t1610);
t1613 = koppelP(4,2);
t1617 = koppelP(4,1);
t1541 = -t1576 * t1617 - t1577 * t1613;
t1545 = -t1576 * t1613 + t1577 * t1617;
t1637 = t1541 * t1572 - t1545 * t1568;
t1614 = koppelP(3,2);
t1618 = koppelP(3,1);
t1542 = -t1576 * t1618 - t1577 * t1614;
t1546 = -t1576 * t1614 + t1577 * t1618;
t1636 = t1542 * t1573 - t1546 * t1569;
t1615 = koppelP(2,2);
t1619 = koppelP(2,1);
t1543 = -t1576 * t1619 - t1577 * t1615;
t1547 = -t1576 * t1615 + t1577 * t1619;
t1635 = t1543 * t1574 - t1547 * t1570;
t1616 = koppelP(1,2);
t1620 = koppelP(1,1);
t1544 = -t1576 * t1620 - t1577 * t1616;
t1548 = -t1576 * t1616 + t1577 * t1620;
t1634 = t1544 * t1575 - t1548 * t1571;
t1633 = pkin(3) * t1683 - t1557 * t1586;
t1632 = pkin(3) * t1680 - t1559 * t1586;
t1631 = pkin(3) * t1678 - t1560 * t1586;
t1630 = pkin(3) * t1676 - t1561 * t1586;
t1612 = rSges(4,1);
t1611 = rSges(4,2);
t1537 = t1583 * t1662 - t1585 * t1606;
t1536 = t1583 * t1664 - t1585 * t1604;
t1535 = t1583 * t1666 - t1585 * t1602;
t1530 = t1583 * t1669 - t1585 * t1590;
t1516 = t1564 * t1585 + t1630 * t1583;
t1515 = t1563 * t1585 + t1631 * t1583;
t1514 = t1562 * t1585 + t1632 * t1583;
t1513 = t1558 * t1585 + t1633 * t1583;
t1504 = -(t1537 * t1575 - t1571 * t1675) * t1720 + (t1516 * t1575 + t1524 * t1571) * t1605 + (t1586 * t1571 + t1575 * t1685) * t1724;
t1503 = (t1537 * t1571 + t1575 * t1675) * t1720 + (-t1516 * t1571 + t1524 * t1575) * t1605 + (-t1571 * t1685 + t1586 * t1575) * t1724;
t1502 = -(t1536 * t1574 - t1570 * t1677) * t1721 + (t1515 * t1574 + t1523 * t1570) * t1603 + (t1586 * t1570 + t1574 * t1685) * t1725;
t1501 = (t1536 * t1570 + t1574 * t1677) * t1721 + (-t1515 * t1570 + t1523 * t1574) * t1603 + (-t1570 * t1685 + t1586 * t1574) * t1725;
t1500 = -(t1535 * t1573 - t1569 * t1679) * t1722 + (t1514 * t1573 + t1522 * t1569) * t1601 + (t1586 * t1569 + t1573 * t1685) * t1726;
t1499 = (t1535 * t1569 + t1573 * t1679) * t1722 + (-t1514 * t1569 + t1522 * t1573) * t1601 + (-t1569 * t1685 + t1586 * t1573) * t1726;
t1498 = -(t1530 * t1572 - t1568 * t1682) * t1723 + (t1513 * t1572 + t1521 * t1568) * t1589 + (t1586 * t1568 + t1572 * t1685) * t1727;
t1497 = (t1530 * t1568 + t1572 * t1682) * t1723 + (-t1513 * t1568 + t1521 * t1572) * t1589 + (-t1568 * t1685 + t1586 * t1572) * t1727;
t1 = [-t1572 * t1653 - t1573 * t1652 - t1574 * t1651 - t1575 * t1650 - m(4) * g(1) + (-t1498 * t1706 - t1500 * t1705 - t1502 * t1704 - t1504 * t1703) * t1579 + (t1572 * t1657 + t1573 * t1656 + t1574 * t1655 + t1575 * t1654) * t1728; t1568 * t1653 + t1569 * t1652 + t1570 * t1651 + t1571 * t1650 - m(4) * g(2) + (-t1497 * t1706 - t1499 * t1705 - t1501 * t1704 - t1503 * t1703) * t1579 + (-t1568 * t1657 - t1569 * t1656 - t1570 * t1655 - t1571 * t1654) * t1728; -m(4) * g(3) + (-(-t1540 * t1720 - t1564 * t1583 * t1605 + (pkin(2) * t1676 + t1630 * t1605) * t1585) * t1693 + (t1599 * t1537 + t1583 * t1672) * t1496) * t1512 + (-(-t1539 * t1721 - t1563 * t1583 * t1603 + (pkin(2) * t1678 + t1631 * t1603) * t1585) * t1696 + (t1597 * t1536 + t1583 * t1673) * t1495) * t1511 + (-(-t1538 * t1722 - t1562 * t1583 * t1601 + (pkin(2) * t1680 + t1632 * t1601) * t1585) * t1699 + (t1595 * t1535 + t1583 * t1674) * t1494) * t1510 + (-(-t1531 * t1723 - t1558 * t1583 * t1589 + (pkin(2) * t1683 + t1633 * t1589) * t1585) * t1702 + (t1587 * t1530 + t1583 * t1681) * t1493) * t1509 + (((t1583 * t1659 + t1585 * t1600) * t1716 + t1564 * t1684 + t1585 * t1561) * t1711 + ((t1583 * t1660 + t1585 * t1598) * t1717 + t1563 * t1684 + t1585 * t1560) * t1712 + ((t1583 * t1661 + t1585 * t1596) * t1718 + t1562 * t1684 + t1585 * t1559) * t1713 + ((t1583 * t1668 + t1585 * t1588) * t1719 + t1558 * t1684 + t1585 * t1557) * t1714) * t1728; ((g(1) * t1612 + g(2) * t1611) * t1576 + (g(1) * t1611 - g(2) * t1612) * t1577) * m(4) + (-(t1503 * t1548 + t1504 * t1544) * t1693 - t1634 * t1707) * t1512 + (-(t1501 * t1547 + t1502 * t1543) * t1696 - t1635 * t1708) * t1511 + (-(t1499 * t1546 + t1500 * t1542) * t1699 - t1636 * t1709) * t1510 + (-(t1497 * t1545 + t1498 * t1541) * t1702 - t1637 * t1710) * t1509 + (t1634 * t1654 + t1635 * t1655 + t1636 * t1656 + t1637 * t1657) * t1728;];
taugX  = t1;
