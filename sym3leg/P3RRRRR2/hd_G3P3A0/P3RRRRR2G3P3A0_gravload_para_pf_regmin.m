% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR2G3P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR2G3P3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:12:33
% EndTime: 2020-03-09 21:12:34
% DurationCPUTime: 1.07s
% Computational Cost: add. (717->147), mult. (1209->310), div. (303->11), fcn. (1404->30), ass. (0->162)
t1650 = cos(qJ(3,1));
t1733 = t1650 ^ 2;
t1647 = cos(qJ(3,2));
t1732 = t1647 ^ 2;
t1644 = cos(qJ(3,3));
t1731 = t1644 ^ 2;
t1632 = legFrame(3,2);
t1611 = sin(t1632);
t1614 = cos(t1632);
t1598 = t1614 * g(1) - t1611 * g(2);
t1629 = qJ(1,3) + qJ(2,3);
t1605 = sin(t1629);
t1608 = cos(t1629);
t1577 = g(3) * t1608 + t1598 * t1605;
t1636 = sin(qJ(2,3));
t1617 = 0.1e1 / t1636;
t1722 = t1577 * t1617;
t1634 = legFrame(1,2);
t1613 = sin(t1634);
t1616 = cos(t1634);
t1600 = t1616 * g(1) - t1613 * g(2);
t1631 = qJ(1,1) + qJ(2,1);
t1607 = sin(t1631);
t1610 = cos(t1631);
t1581 = g(3) * t1610 + t1600 * t1607;
t1642 = sin(qJ(2,1));
t1619 = 0.1e1 / t1642;
t1721 = t1581 * t1619;
t1633 = legFrame(2,2);
t1612 = sin(t1633);
t1615 = cos(t1633);
t1599 = t1615 * g(1) - t1612 * g(2);
t1630 = qJ(1,2) + qJ(2,2);
t1606 = sin(t1630);
t1609 = cos(t1630);
t1661 = g(3) * t1609 + t1599 * t1606;
t1639 = sin(qJ(2,2));
t1618 = 0.1e1 / t1639;
t1624 = 0.1e1 / t1647;
t1698 = t1618 * t1624;
t1730 = t1661 * t1698;
t1714 = t1606 * t1618;
t1729 = t1661 * t1714;
t1646 = cos(qJ(1,3));
t1728 = pkin(1) * t1646;
t1649 = cos(qJ(1,2));
t1727 = pkin(1) * t1649;
t1652 = cos(qJ(1,1));
t1726 = pkin(1) * t1652;
t1723 = t1661 * t1618;
t1637 = sin(qJ(1,3));
t1645 = cos(qJ(2,3));
t1592 = t1637 * t1636 - t1646 * t1645;
t1720 = t1592 * t1644;
t1640 = sin(qJ(1,2));
t1648 = cos(qJ(2,2));
t1593 = t1640 * t1639 - t1649 * t1648;
t1719 = t1593 * t1647;
t1643 = sin(qJ(1,1));
t1651 = cos(qJ(2,1));
t1594 = t1643 * t1642 - t1652 * t1651;
t1718 = t1594 * t1650;
t1715 = t1605 * t1617;
t1713 = t1607 * t1619;
t1621 = 0.1e1 / t1644;
t1712 = t1611 * t1621;
t1635 = sin(qJ(3,3));
t1711 = t1611 * t1635;
t1710 = t1612 * t1624;
t1638 = sin(qJ(3,2));
t1709 = t1612 * t1638;
t1627 = 0.1e1 / t1650;
t1708 = t1613 * t1627;
t1641 = sin(qJ(3,1));
t1707 = t1613 * t1641;
t1706 = t1614 * t1621;
t1705 = t1614 * t1635;
t1704 = t1615 * t1624;
t1703 = t1615 * t1638;
t1702 = t1616 * t1627;
t1701 = t1616 * t1641;
t1700 = t1617 * t1621;
t1622 = 0.1e1 / t1731;
t1699 = t1617 * t1622;
t1625 = 0.1e1 / t1732;
t1697 = t1618 * t1625;
t1696 = t1619 * t1627;
t1628 = 0.1e1 / t1733;
t1695 = t1619 * t1628;
t1694 = pkin(2) * t1592 * t1731;
t1693 = pkin(2) * t1593 * t1732;
t1692 = pkin(2) * t1594 * t1733;
t1691 = t1645 * t1635 * pkin(1);
t1690 = t1648 * t1638 * pkin(1);
t1689 = t1651 * t1641 * pkin(1);
t1687 = t1638 * t1723;
t1686 = t1635 * t1722;
t1685 = t1641 * t1721;
t1684 = t1577 * t1715;
t1683 = t1577 * t1700;
t1682 = t1577 * t1699;
t1578 = -g(3) * t1605 + t1598 * t1608;
t1681 = t1578 * t1700;
t1680 = t1578 * t1699;
t1678 = t1661 * t1697;
t1580 = -g(3) * t1606 + t1599 * t1609;
t1677 = t1580 * t1698;
t1676 = t1580 * t1697;
t1675 = t1581 * t1713;
t1674 = t1581 * t1696;
t1673 = t1581 * t1695;
t1582 = -g(3) * t1607 + t1600 * t1610;
t1672 = t1582 * t1696;
t1671 = t1582 * t1695;
t1583 = pkin(2) * (t1646 * t1636 + t1637 * t1645) * t1644 + t1637 * pkin(1);
t1670 = t1583 * t1700;
t1584 = pkin(2) * (t1649 * t1639 + t1640 * t1648) * t1647 + t1640 * pkin(1);
t1669 = t1584 * t1698;
t1585 = pkin(2) * (t1652 * t1642 + t1643 * t1651) * t1650 + t1643 * pkin(1);
t1668 = t1585 * t1696;
t1586 = g(3) * t1646 + t1598 * t1637;
t1667 = t1586 * t1700;
t1587 = -g(3) * t1637 + t1598 * t1646;
t1666 = t1587 * t1700;
t1588 = g(3) * t1649 + t1599 * t1640;
t1665 = t1588 * t1698;
t1589 = -g(3) * t1640 + t1599 * t1649;
t1664 = t1589 * t1698;
t1590 = g(3) * t1652 + t1600 * t1643;
t1663 = t1590 * t1696;
t1591 = -g(3) * t1643 + t1600 * t1652;
t1662 = t1591 * t1696;
t1660 = t1624 * t1687;
t1659 = t1625 * t1687;
t1658 = t1621 * t1686;
t1657 = t1622 * t1686;
t1656 = t1627 * t1685;
t1655 = t1628 * t1685;
t1654 = 0.1e1 / pkin(1);
t1653 = 0.1e1 / pkin(2);
t1597 = t1613 * g(1) + t1616 * g(2);
t1596 = t1612 * g(1) + t1615 * g(2);
t1595 = t1611 * g(1) + t1614 * g(2);
t1573 = -t1616 * t1718 + t1707;
t1572 = t1613 * t1718 + t1701;
t1571 = -t1615 * t1719 + t1709;
t1570 = t1612 * t1719 + t1703;
t1569 = -t1614 * t1720 + t1711;
t1568 = t1611 * t1720 + t1705;
t1567 = t1582 * t1650 + t1597 * t1641;
t1566 = t1582 * t1641 - t1597 * t1650;
t1565 = t1580 * t1647 + t1596 * t1638;
t1564 = t1580 * t1638 - t1596 * t1647;
t1563 = t1578 * t1644 + t1595 * t1635;
t1562 = t1578 * t1635 - t1595 * t1644;
t1561 = t1616 * t1692 + (-pkin(2) * t1707 - t1616 * t1726) * t1650 - t1613 * t1689;
t1560 = t1615 * t1693 + (-pkin(2) * t1709 - t1615 * t1727) * t1647 - t1612 * t1690;
t1559 = t1614 * t1694 + (-pkin(2) * t1711 - t1614 * t1728) * t1644 - t1611 * t1691;
t1558 = -t1613 * t1692 + (-pkin(2) * t1701 + t1613 * t1726) * t1650 - t1616 * t1689;
t1557 = -t1612 * t1693 + (-pkin(2) * t1703 + t1612 * t1727) * t1647 - t1615 * t1690;
t1556 = -t1611 * t1694 + (-pkin(2) * t1705 + t1611 * t1728) * t1644 - t1614 * t1691;
t1 = [0, (t1569 * t1667 + t1571 * t1665 + t1573 * t1663) * t1654, (t1569 * t1666 + t1571 * t1664 + t1573 * t1662) * t1654, 0, (t1569 * t1683 + t1571 * t1730 + t1573 * t1674 + (t1559 * t1682 + t1560 * t1678 + t1561 * t1673) * t1653) * t1654, (t1569 * t1681 + t1571 * t1677 + t1573 * t1672 + (t1559 * t1680 + t1560 * t1676 + t1561 * t1671) * t1653) * t1654, 0, 0, 0, 0, 0, (t1569 * t1722 + t1571 * t1723 + t1573 * t1721) * t1654 + (t1562 * t1712 + t1564 * t1710 + t1566 * t1708 + (t1559 * t1683 + t1560 * t1730 + t1561 * t1674) * t1654) * t1653, (-t1569 * t1658 - t1571 * t1660 - t1573 * t1656) * t1654 + (t1563 * t1712 + t1565 * t1710 + t1567 * t1708 + (-t1559 * t1657 - t1560 * t1659 - t1561 * t1655) * t1654) * t1653, -g(1); 0, (t1568 * t1667 + t1570 * t1665 + t1572 * t1663) * t1654, (t1568 * t1666 + t1570 * t1664 + t1572 * t1662) * t1654, 0, (t1568 * t1683 + t1570 * t1730 + t1572 * t1674 + (t1556 * t1682 + t1557 * t1678 + t1558 * t1673) * t1653) * t1654, (t1568 * t1681 + t1570 * t1677 + t1572 * t1672 + (t1556 * t1680 + t1557 * t1676 + t1558 * t1671) * t1653) * t1654, 0, 0, 0, 0, 0, (t1568 * t1722 + t1570 * t1723 + t1572 * t1721) * t1654 + (t1562 * t1706 + t1564 * t1704 + t1566 * t1702 + (t1556 * t1683 + t1557 * t1730 + t1558 * t1674) * t1654) * t1653, (-t1568 * t1658 - t1570 * t1660 - t1572 * t1656) * t1654 + (t1563 * t1706 + t1565 * t1704 + t1567 * t1702 + (-t1556 * t1657 - t1557 * t1659 - t1558 * t1655) * t1654) * t1653, -g(2); 0, (-t1586 * t1715 - t1588 * t1714 - t1590 * t1713) * t1654, (-t1587 * t1715 - t1589 * t1714 - t1591 * t1713) * t1654, 0, (-t1684 - t1729 - t1675 + (t1577 * t1670 + t1581 * t1668 + t1661 * t1669) * t1653) * t1654, (-t1578 * t1715 - t1580 * t1714 - t1582 * t1713 + (t1578 * t1670 + t1580 * t1669 + t1582 * t1668) * t1653) * t1654, 0, 0, 0, 0, 0, (-t1647 * t1729 - t1644 * t1684 - t1650 * t1675 + (t1583 * t1722 + t1584 * t1723 + t1585 * t1721) * t1653) * t1654, (t1606 * t1687 + t1605 * t1686 + t1607 * t1685 + (-t1583 * t1658 - t1584 * t1660 - t1585 * t1656) * t1653) * t1654, -g(3);];
tau_reg  = t1;
