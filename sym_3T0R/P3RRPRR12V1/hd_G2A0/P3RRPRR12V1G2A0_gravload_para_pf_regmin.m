% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V1G2A0
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:06:24
% EndTime: 2020-08-06 19:06:26
% DurationCPUTime: 1.44s
% Computational Cost: add. (1035->148), mult. (1881->282), div. (171->6), fcn. (1761->18), ass. (0->152)
t1652 = legFrame(3,2);
t1640 = sin(t1652);
t1643 = cos(t1652);
t1613 = t1643 * g(1) - t1640 * g(2);
t1656 = sin(qJ(1,3));
t1662 = cos(qJ(1,3));
t1604 = -g(3) * t1656 + t1613 * t1662;
t1726 = t1604 * t1662;
t1653 = legFrame(2,2);
t1641 = sin(t1653);
t1644 = cos(t1653);
t1614 = t1644 * g(1) - t1641 * g(2);
t1658 = sin(qJ(1,2));
t1664 = cos(qJ(1,2));
t1605 = -g(3) * t1658 + t1614 * t1664;
t1724 = t1605 * t1664;
t1654 = legFrame(1,2);
t1642 = sin(t1654);
t1645 = cos(t1654);
t1615 = t1645 * g(1) - t1642 * g(2);
t1660 = sin(qJ(1,1));
t1666 = cos(qJ(1,1));
t1606 = -g(3) * t1660 + t1615 * t1666;
t1722 = t1606 * t1666;
t1655 = sin(qJ(2,3));
t1712 = t1655 * t1656;
t1754 = t1662 * pkin(4);
t1616 = qJ(3,3) * t1712 + t1754;
t1661 = cos(qJ(2,3));
t1649 = t1661 ^ 2;
t1667 = pkin(1) + pkin(2);
t1703 = t1667 * t1655;
t1710 = t1656 * t1667;
t1748 = t1643 * qJ(3,3);
t1577 = (-t1640 * t1710 - t1748) * t1649 + (-t1640 * t1616 + t1643 * t1703) * t1661 + t1748;
t1657 = sin(qJ(2,2));
t1709 = t1657 * t1658;
t1753 = t1664 * pkin(4);
t1617 = qJ(3,2) * t1709 + t1753;
t1663 = cos(qJ(2,2));
t1650 = t1663 ^ 2;
t1702 = t1667 * t1657;
t1707 = t1658 * t1667;
t1747 = t1644 * qJ(3,2);
t1578 = (-t1641 * t1707 - t1747) * t1650 + (-t1641 * t1617 + t1644 * t1702) * t1663 + t1747;
t1659 = sin(qJ(2,1));
t1706 = t1659 * t1660;
t1752 = t1666 * pkin(4);
t1618 = qJ(3,1) * t1706 + t1752;
t1665 = cos(qJ(2,1));
t1651 = t1665 ^ 2;
t1701 = t1667 * t1659;
t1704 = t1660 * t1667;
t1746 = t1645 * qJ(3,1);
t1579 = (-t1642 * t1704 - t1746) * t1651 + (-t1642 * t1618 + t1645 * t1701) * t1665 + t1746;
t1637 = t1655 * qJ(3,3);
t1769 = t1667 * t1661 + t1637;
t1619 = 0.1e1 / t1769;
t1638 = t1657 * qJ(3,2);
t1768 = t1667 * t1663 + t1638;
t1620 = 0.1e1 / t1768;
t1639 = t1659 * qJ(3,1);
t1767 = t1667 * t1665 + t1639;
t1621 = 0.1e1 / t1767;
t1612 = t1642 * g(1) + t1645 * g(2);
t1784 = g(3) * t1666 + t1615 * t1660;
t1589 = t1612 * t1659 + t1665 * t1784;
t1670 = 0.1e1 / qJ(3,1);
t1791 = t1589 * t1670;
t1611 = t1641 * g(1) + t1644 * g(2);
t1783 = g(3) * t1664 + t1614 * t1658;
t1585 = t1611 * t1657 + t1663 * t1783;
t1669 = 0.1e1 / qJ(3,2);
t1792 = t1585 * t1669;
t1610 = t1640 * g(1) + t1643 * g(2);
t1782 = g(3) * t1662 + t1613 * t1656;
t1581 = t1610 * t1655 + t1661 * t1782;
t1668 = 0.1e1 / qJ(3,3);
t1793 = t1581 * t1668;
t1801 = t1659 * t1722;
t1802 = t1657 * t1724;
t1803 = t1655 * t1726;
t1805 = (t1577 * t1793 - t1640 * t1803) * t1619 + (t1578 * t1792 - t1641 * t1802) * t1620 + (t1579 * t1791 - t1642 * t1801) * t1621;
t1751 = t1640 * qJ(3,3);
t1574 = (t1643 * t1710 - t1751) * t1649 + (t1616 * t1643 + t1640 * t1703) * t1661 + t1751;
t1750 = t1641 * qJ(3,2);
t1575 = (t1644 * t1707 - t1750) * t1650 + (t1617 * t1644 + t1641 * t1702) * t1663 + t1750;
t1749 = t1642 * qJ(3,1);
t1576 = (t1645 * t1704 - t1749) * t1651 + (t1618 * t1645 + t1642 * t1701) * t1665 + t1749;
t1804 = (t1574 * t1793 + t1643 * t1803) * t1619 + (t1575 * t1792 + t1644 * t1802) * t1620 + (t1576 * t1791 + t1645 * t1801) * t1621;
t1779 = -t1660 * pkin(4) + t1666 * t1767;
t1719 = t1779 * t1670;
t1689 = t1665 * t1719;
t1780 = -t1658 * pkin(4) + t1664 * t1768;
t1720 = t1780 * t1669;
t1690 = t1663 * t1720;
t1781 = -t1656 * pkin(4) + t1662 * t1769;
t1721 = t1781 * t1668;
t1691 = t1661 * t1721;
t1794 = (t1589 * t1689 - t1606 * t1706) * t1621 + (t1585 * t1690 - t1605 * t1709) * t1620 + (t1581 * t1691 - t1604 * t1712) * t1619;
t1580 = -t1610 * t1661 + t1655 * t1782;
t1775 = t1580 * t1668;
t1584 = -t1611 * t1663 + t1657 * t1783;
t1774 = t1584 * t1669;
t1588 = -t1612 * t1665 + t1659 * t1784;
t1773 = t1588 * t1670;
t1745 = t1661 * qJ(3,3);
t1744 = t1663 * qJ(3,2);
t1743 = t1665 * qJ(3,1);
t1631 = t1661 * pkin(1) + t1637;
t1571 = -t1610 * t1631 + t1782 * (t1655 * pkin(1) - t1745);
t1742 = t1571 * t1668;
t1632 = t1663 * pkin(1) + t1638;
t1572 = -t1611 * t1632 + t1783 * (t1657 * pkin(1) - t1744);
t1741 = t1572 * t1669;
t1633 = t1665 * pkin(1) + t1639;
t1573 = -t1612 * t1633 + t1784 * (t1659 * pkin(1) - t1743);
t1740 = t1573 * t1670;
t1727 = t1604 * t1656;
t1725 = t1605 * t1658;
t1723 = t1606 * t1660;
t1718 = t1619 * t1656;
t1717 = t1619 * t1662;
t1716 = t1620 * t1658;
t1715 = t1620 * t1664;
t1714 = t1621 * t1660;
t1713 = t1621 * t1666;
t1697 = t1631 * t1726;
t1696 = t1661 * t1726;
t1695 = t1632 * t1724;
t1694 = t1663 * t1724;
t1693 = t1633 * t1722;
t1692 = t1665 * t1722;
t1688 = t1640 * t1717;
t1687 = t1643 * t1717;
t1686 = t1641 * t1715;
t1685 = t1644 * t1715;
t1684 = t1642 * t1713;
t1683 = t1645 * t1713;
t1673 = t1714 * t1784 + t1716 * t1783 + t1718 * t1782;
t1672 = t1684 * t1784 + t1686 * t1783 + t1688 * t1782;
t1671 = t1683 * t1784 + t1685 * t1783 + t1687 * t1782;
t1624 = t1701 - t1743;
t1623 = t1702 - t1744;
t1622 = t1703 - t1745;
t1600 = t1660 * t1767 + t1752;
t1599 = t1658 * t1768 + t1753;
t1598 = t1656 * t1769 + t1754;
t1570 = (t1588 * t1719 + t1723) * t1665 * t1621 + (t1584 * t1720 + t1725) * t1663 * t1620 + (t1580 * t1721 + t1727) * t1661 * t1619;
t1569 = (t1576 * t1773 - t1645 * t1692) * t1621 + (t1575 * t1774 - t1644 * t1694) * t1620 + (t1574 * t1775 - t1643 * t1696) * t1619;
t1568 = (t1579 * t1773 + t1642 * t1692) * t1621 + (t1578 * t1774 + t1641 * t1694) * t1620 + (t1577 * t1775 + t1640 * t1696) * t1619;
t1 = [0, -t1604 * t1687 - t1605 * t1685 - t1606 * t1683, t1671, 0, 0, 0, 0, 0, t1569, t1804, t1569, -t1671, -t1804, -(t1600 * t1645 + t1642 * t1624) * t1773 - (t1599 * t1644 + t1641 * t1623) * t1774 - (t1598 * t1643 + t1640 * t1622) * t1775 + (t1576 * t1740 - t1645 * t1693) * t1621 + (t1575 * t1741 - t1644 * t1695) * t1620 + (t1574 * t1742 - t1643 * t1697) * t1619, -g(1); 0, t1604 * t1688 + t1605 * t1686 + t1606 * t1684, -t1672, 0, 0, 0, 0, 0, t1568, t1805, t1568, t1672, -t1805, -(-t1600 * t1642 + t1624 * t1645) * t1773 - (-t1599 * t1641 + t1623 * t1644) * t1774 - (-t1598 * t1640 + t1622 * t1643) * t1775 + (t1579 * t1740 + t1642 * t1693) * t1621 + (t1578 * t1741 + t1641 * t1695) * t1620 + (t1577 * t1742 + t1640 * t1697) * t1619, -g(2); 0, t1604 * t1718 + t1605 * t1716 + t1606 * t1714, -t1673, 0, 0, 0, 0, 0, t1570, t1794, t1570, t1673, -t1794, -t1779 * t1773 - t1780 * t1774 - t1781 * t1775 + (t1573 * t1689 + t1633 * t1723) * t1621 + (t1572 * t1690 + t1632 * t1725) * t1620 + (t1571 * t1691 + t1631 * t1727) * t1619, -g(3);];
tau_reg  = t1;
