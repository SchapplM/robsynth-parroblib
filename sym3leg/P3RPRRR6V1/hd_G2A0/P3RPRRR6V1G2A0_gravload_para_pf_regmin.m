% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR6V1G2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR6V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:37:13
% EndTime: 2020-08-06 18:37:14
% DurationCPUTime: 1.43s
% Computational Cost: add. (900->187), mult. (1059->292), div. (111->10), fcn. (948->53), ass. (0->164)
t1705 = sin(pkin(7));
t1723 = -pkin(6) - pkin(5);
t1665 = t1723 * t1705 - pkin(1);
t1711 = sin(qJ(1,3));
t1656 = t1665 * t1711;
t1706 = cos(pkin(7));
t1717 = cos(qJ(1,3));
t1765 = t1717 * t1723;
t1766 = t1717 * t1705;
t1817 = t1656 - pkin(2) * t1766 - (t1711 * pkin(2) + t1765) * t1706;
t1713 = sin(qJ(1,2));
t1657 = t1665 * t1713;
t1719 = cos(qJ(1,2));
t1763 = t1719 * t1723;
t1764 = t1719 * t1705;
t1816 = t1657 - pkin(2) * t1764 - (t1713 * pkin(2) + t1763) * t1706;
t1715 = sin(qJ(1,1));
t1658 = t1665 * t1715;
t1721 = cos(qJ(1,1));
t1761 = t1721 * t1723;
t1762 = t1721 * t1705;
t1815 = t1658 - pkin(2) * t1762 - (t1715 * pkin(2) + t1761) * t1706;
t1707 = legFrame(3,2);
t1684 = sin(t1707);
t1687 = cos(t1707);
t1653 = t1687 * g(1) - t1684 * g(2);
t1638 = g(3) * t1717 + t1653 * t1711;
t1814 = t1638 / 0.2e1;
t1708 = legFrame(2,2);
t1685 = sin(t1708);
t1688 = cos(t1708);
t1654 = t1688 * g(1) - t1685 * g(2);
t1639 = g(3) * t1719 + t1654 * t1713;
t1813 = t1639 / 0.2e1;
t1709 = legFrame(1,2);
t1686 = sin(t1709);
t1689 = cos(t1709);
t1655 = t1689 * g(1) - t1686 * g(2);
t1640 = g(3) * t1721 + t1655 * t1715;
t1812 = t1640 / 0.2e1;
t1720 = cos(qJ(3,1));
t1676 = t1720 * pkin(3) + pkin(2);
t1683 = t1706 * pkin(1);
t1661 = 0.1e1 / (t1683 + t1676);
t1698 = qJ(1,1) + pkin(7);
t1679 = sin(t1698);
t1682 = cos(t1698);
t1788 = (-g(3) * t1679 + t1655 * t1682) * t1661;
t1743 = t1720 * t1788;
t1811 = -t1743 / 0.2e1;
t1714 = sin(qJ(3,1));
t1744 = t1714 * t1788;
t1810 = t1744 / 0.2e1;
t1718 = cos(qJ(3,2));
t1675 = t1718 * pkin(3) + pkin(2);
t1660 = 0.1e1 / (t1683 + t1675);
t1697 = qJ(1,2) + pkin(7);
t1678 = sin(t1697);
t1681 = cos(t1697);
t1789 = (-g(3) * t1678 + t1654 * t1681) * t1660;
t1745 = t1718 * t1789;
t1809 = -t1745 / 0.2e1;
t1712 = sin(qJ(3,2));
t1746 = t1712 * t1789;
t1808 = t1746 / 0.2e1;
t1716 = cos(qJ(3,3));
t1674 = t1716 * pkin(3) + pkin(2);
t1659 = 0.1e1 / (t1683 + t1674);
t1696 = qJ(1,3) + pkin(7);
t1677 = sin(t1696);
t1680 = cos(t1696);
t1790 = (-g(3) * t1677 + t1653 * t1680) * t1659;
t1747 = t1716 * t1790;
t1807 = -t1747 / 0.2e1;
t1710 = sin(qJ(3,3));
t1748 = t1710 * t1790;
t1806 = t1748 / 0.2e1;
t1737 = g(3) * t1715 - t1655 * t1721;
t1671 = t1709 + t1698;
t1672 = -t1709 + t1698;
t1646 = cos(t1672) + cos(t1671);
t1782 = t1646 * t1661;
t1805 = t1737 * t1782;
t1643 = -sin(t1671) + sin(t1672);
t1785 = t1643 * t1661;
t1804 = t1737 * t1785;
t1738 = g(3) * t1713 - t1654 * t1719;
t1669 = t1708 + t1697;
t1670 = -t1708 + t1697;
t1645 = cos(t1670) + cos(t1669);
t1783 = t1645 * t1660;
t1803 = t1738 * t1783;
t1642 = -sin(t1669) + sin(t1670);
t1786 = t1642 * t1660;
t1802 = t1738 * t1786;
t1739 = g(3) * t1711 - t1653 * t1717;
t1667 = t1707 + t1696;
t1668 = -t1707 + t1696;
t1644 = cos(t1668) + cos(t1667);
t1784 = t1644 * t1659;
t1801 = t1739 * t1784;
t1641 = -sin(t1667) + sin(t1668);
t1787 = t1641 * t1659;
t1800 = t1739 * t1787;
t1777 = t1661 * t1679;
t1779 = t1660 * t1678;
t1781 = t1659 * t1677;
t1799 = -t1737 * t1777 - t1738 * t1779 - t1739 * t1781;
t1798 = -0.2e1 * pkin(1);
t1797 = -0.2e1 * pkin(2);
t1796 = 0.2e1 * pkin(2);
t1795 = 0.2e1 * t1723;
t1794 = pkin(1) / 0.2e1;
t1755 = pkin(7) + qJ(3,3);
t1758 = -pkin(7) + qJ(3,3);
t1793 = (t1677 * t1795 + t1717 * t1798 + t1680 * t1797 + (-cos(qJ(1,3) - t1758) - cos(qJ(1,3) + t1755)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t1710 * t1796 + (sin(t1755) + sin(t1758)) * pkin(1));
t1756 = pkin(7) + qJ(3,2);
t1759 = -pkin(7) + qJ(3,2);
t1792 = (t1678 * t1795 + t1719 * t1798 + t1681 * t1797 + (-cos(qJ(1,2) - t1759) - cos(qJ(1,2) + t1756)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t1712 * t1796 + (sin(t1756) + sin(t1759)) * pkin(1));
t1757 = pkin(7) + qJ(3,1);
t1760 = -pkin(7) + qJ(3,1);
t1791 = (t1679 * t1795 + t1721 * t1798 + t1682 * t1797 + (-cos(qJ(1,1) - t1760) - cos(qJ(1,1) + t1757)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t1714 * t1796 + (sin(t1757) + sin(t1760)) * pkin(1));
t1780 = t1659 / t1710;
t1778 = t1660 / t1712;
t1776 = t1661 / t1714;
t1775 = t1674 * t1711;
t1774 = t1675 * t1713;
t1773 = t1676 * t1715;
t1772 = t1710 * t1684;
t1771 = t1710 * t1687;
t1770 = t1712 * t1685;
t1769 = t1712 * t1688;
t1768 = t1714 * t1686;
t1767 = t1714 * t1689;
t1754 = (t1711 * t1706 + t1766) * t1716 ^ 2 * pkin(3);
t1753 = (t1713 * t1706 + t1764) * t1718 ^ 2 * pkin(3);
t1752 = (t1715 * t1706 + t1762) * t1720 ^ 2 * pkin(3);
t1751 = ((t1765 + t1775) * t1706 - t1656 + t1674 * t1766) * t1780;
t1750 = ((t1763 + t1774) * t1706 - t1657 + t1675 * t1764) * t1778;
t1749 = ((t1761 + t1773) * t1706 - t1658 + t1676 * t1762) * t1776;
t1650 = t1684 * g(1) + t1687 * g(2);
t1742 = t1650 * t1780;
t1651 = t1685 * g(1) + t1688 * g(2);
t1741 = t1651 * t1778;
t1652 = t1686 * g(1) + t1689 * g(2);
t1740 = t1652 * t1776;
t1736 = t1684 * t1751;
t1735 = t1687 * t1751;
t1734 = t1685 * t1750;
t1733 = t1688 * t1750;
t1732 = t1686 * t1749;
t1731 = t1689 * t1749;
t1730 = g(3) * t1680 + t1653 * t1677;
t1729 = g(3) * t1681 + t1654 * t1678;
t1728 = g(3) * t1682 + t1655 * t1679;
t1724 = 0.1e1 / pkin(3);
t1673 = t1683 + pkin(2);
t1619 = t1652 * t1714 + t1728 * t1720;
t1618 = -t1652 * t1720 + t1728 * t1714;
t1617 = t1651 * t1712 + t1729 * t1718;
t1616 = -t1651 * t1718 + t1729 * t1712;
t1615 = t1650 * t1710 + t1730 * t1716;
t1614 = -t1650 * t1716 + t1730 * t1710;
t1 = [0, t1805 / 0.2e1 + t1803 / 0.2e1 + t1801 / 0.2e1, t1782 * t1812 + t1783 * t1813 + t1784 * t1814, -(t1689 * t1752 + (pkin(3) * t1768 - t1815 * t1689) * t1720 + t1673 * t1768) * t1740 - (t1688 * t1753 + (pkin(3) * t1770 - t1816 * t1688) * t1718 + t1673 * t1770) * t1741 - (t1687 * t1754 + (pkin(3) * t1772 - t1817 * t1687) * t1716 + t1673 * t1772) * t1742 + (t1801 + t1803 + t1805) * t1794, 0, 0, 0, 0, 0, (-t1614 * t1735 - t1616 * t1733 - t1618 * t1731) * t1724 + t1644 * t1807 + t1645 * t1809 + t1646 * t1811, (-t1615 * t1735 - t1617 * t1733 - t1619 * t1731) * t1724 + t1644 * t1806 + t1645 * t1808 + t1646 * t1810, -g(1); 0, t1804 / 0.2e1 + t1802 / 0.2e1 + t1800 / 0.2e1, t1785 * t1812 + t1786 * t1813 + t1787 * t1814, -(-t1686 * t1752 + (pkin(3) * t1767 + t1815 * t1686) * t1720 + t1673 * t1767) * t1740 - (-t1685 * t1753 + (pkin(3) * t1769 + t1816 * t1685) * t1718 + t1673 * t1769) * t1741 - (-t1684 * t1754 + (pkin(3) * t1771 + t1817 * t1684) * t1716 + t1673 * t1771) * t1742 + (t1800 + t1802 + t1804) * t1794, 0, 0, 0, 0, 0, (t1614 * t1736 + t1616 * t1734 + t1618 * t1732) * t1724 + t1641 * t1807 + t1642 * t1809 + t1643 * t1811, (t1615 * t1736 + t1617 * t1734 + t1619 * t1732) * t1724 + t1641 * t1806 + t1642 * t1808 + t1643 * t1810, -g(2); 0, t1799, -t1638 * t1781 - t1639 * t1779 - t1640 * t1777, -((t1676 * t1721 - t1715 * t1723) * t1706 - t1665 * t1721 - t1705 * t1773) * t1720 * t1740 - ((t1675 * t1719 - t1713 * t1723) * t1706 - t1665 * t1719 - t1705 * t1774) * t1718 * t1741 - ((t1674 * t1717 - t1711 * t1723) * t1706 - t1665 * t1717 - t1705 * t1775) * t1716 * t1742 + t1799 * pkin(1), 0, 0, 0, 0, 0, t1677 * t1747 + t1678 * t1745 + t1679 * t1743 + (t1614 * t1793 + t1616 * t1792 + t1618 * t1791) * t1724, -t1677 * t1748 - t1678 * t1746 - t1679 * t1744 + (t1615 * t1793 + t1617 * t1792 + t1619 * t1791) * t1724, -g(3);];
tau_reg  = t1;
