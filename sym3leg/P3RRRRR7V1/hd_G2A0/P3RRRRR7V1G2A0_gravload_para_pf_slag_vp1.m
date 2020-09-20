% Calculate Gravitation load for parallel robot
% P3RRRRR7V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:41:50
% EndTime: 2020-08-07 03:41:52
% DurationCPUTime: 1.53s
% Computational Cost: add. (864->236), mult. (1662->330), div. (60->14), fcn. (975->60), ass. (0->152)
t1809 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - m(3) * (rSges(3,3) + pkin(4));
t1771 = sin(qJ(3,1));
t1780 = cos(qJ(3,1));
t1823 = -m(2) * rSges(2,1) - pkin(1) * m(3);
t1869 = m(3) * rSges(3,2);
t1870 = m(3) * rSges(3,1);
t1707 = t1771 * t1869 - t1780 * t1870 + t1823;
t1784 = m(2) * rSges(2,2);
t1713 = t1784 + (rSges(3,1) * t1771 + rSges(3,2) * t1780) * m(3);
t1772 = sin(qJ(2,1));
t1781 = cos(qJ(2,1));
t1876 = -t1707 * t1781 - t1713 * t1772;
t1768 = sin(qJ(3,2));
t1777 = cos(qJ(3,2));
t1706 = t1768 * t1869 - t1777 * t1870 + t1823;
t1712 = t1784 + (rSges(3,1) * t1768 + rSges(3,2) * t1777) * m(3);
t1769 = sin(qJ(2,2));
t1778 = cos(qJ(2,2));
t1875 = -t1706 * t1778 - t1712 * t1769;
t1765 = sin(qJ(3,3));
t1774 = cos(qJ(3,3));
t1705 = t1765 * t1869 - t1774 * t1870 + t1823;
t1711 = t1784 + (rSges(3,1) * t1765 + rSges(3,2) * t1774) * m(3);
t1766 = sin(qJ(2,3));
t1775 = cos(qJ(2,3));
t1874 = -t1705 * t1775 - t1711 * t1766;
t1873 = 2 * pkin(2);
t1786 = (-pkin(5) - pkin(4));
t1872 = 2 * t1786;
t1871 = m(1) * rSges(1,1);
t1866 = m(3) / pkin(2);
t1865 = pkin(2) * t1774;
t1864 = pkin(2) * t1777;
t1863 = pkin(2) * t1780;
t1862 = -qJ(3,1) + qJ(1,1);
t1861 = qJ(3,1) + qJ(1,1);
t1860 = -qJ(3,2) + qJ(1,2);
t1859 = qJ(3,2) + qJ(1,2);
t1858 = -qJ(3,3) + qJ(1,3);
t1857 = qJ(3,3) + qJ(1,3);
t1856 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1855 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1854 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1853 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1852 = qJ(1,3) + 0.2e1 * qJ(2,3);
t1851 = -0.2e1 * qJ(2,3) + qJ(1,3);
t1762 = legFrame(3,2);
t1743 = sin(t1762);
t1746 = cos(t1762);
t1717 = g(1) * t1746 - g(2) * t1743;
t1759 = qJ(2,3) + qJ(3,3);
t1740 = cos(t1759);
t1749 = g(3) * t1871;
t1767 = sin(qJ(1,3));
t1776 = cos(qJ(1,3));
t1808 = t1809 * g(3);
t1850 = (t1808 * t1776 + t1767 * (t1874 * g(3) + t1749) + ((-t1871 - t1874) * t1776 + t1767 * t1809) * t1717) / (pkin(1) * t1775 + pkin(2) * t1740);
t1763 = legFrame(2,2);
t1744 = sin(t1763);
t1747 = cos(t1763);
t1718 = g(1) * t1747 - g(2) * t1744;
t1760 = qJ(2,2) + qJ(3,2);
t1741 = cos(t1760);
t1770 = sin(qJ(1,2));
t1779 = cos(qJ(1,2));
t1849 = (t1808 * t1779 + t1770 * (t1875 * g(3) + t1749) + ((-t1871 - t1875) * t1779 + t1770 * t1809) * t1718) / (pkin(1) * t1778 + pkin(2) * t1741);
t1764 = legFrame(1,2);
t1745 = sin(t1764);
t1748 = cos(t1764);
t1719 = g(1) * t1748 - g(2) * t1745;
t1761 = qJ(2,1) + qJ(3,1);
t1742 = cos(t1761);
t1773 = sin(qJ(1,1));
t1782 = cos(qJ(1,1));
t1848 = (t1808 * t1782 + t1773 * (t1876 * g(3) + t1749) + ((-t1871 - t1876) * t1782 + t1773 * t1809) * t1719) / (pkin(1) * t1781 + pkin(2) * t1742);
t1714 = g(1) * t1743 + g(2) * t1746;
t1737 = sin(t1759);
t1750 = 0.1e1 / t1765;
t1807 = g(3) * t1776 + t1717 * t1767;
t1847 = ((-rSges(3,1) * t1714 + t1807 * rSges(3,2)) * t1740 + t1737 * (t1807 * rSges(3,1) + rSges(3,2) * t1714)) * t1750;
t1715 = g(1) * t1744 + g(2) * t1747;
t1738 = sin(t1760);
t1751 = 0.1e1 / t1768;
t1806 = g(3) * t1779 + t1718 * t1770;
t1846 = ((-rSges(3,1) * t1715 + t1806 * rSges(3,2)) * t1741 + t1738 * (t1806 * rSges(3,1) + rSges(3,2) * t1715)) * t1751;
t1716 = g(1) * t1745 + g(2) * t1748;
t1739 = sin(t1761);
t1752 = 0.1e1 / t1771;
t1805 = g(3) * t1782 + t1719 * t1773;
t1845 = ((-rSges(3,1) * t1716 + t1805 * rSges(3,2)) * t1742 + t1739 * (t1805 * rSges(3,1) + rSges(3,2) * t1716)) * t1752;
t1753 = t1774 ^ 2;
t1723 = pkin(1) * t1774 + t1753 * t1873 - pkin(2);
t1838 = t1723 * t1767;
t1755 = t1777 ^ 2;
t1724 = pkin(1) * t1777 + t1755 * t1873 - pkin(2);
t1837 = t1724 * t1770;
t1757 = t1780 ^ 2;
t1725 = pkin(1) * t1780 + t1757 * t1873 - pkin(2);
t1836 = t1725 * t1773;
t1835 = (pkin(1) + 0.2e1 * t1865) * t1765;
t1834 = (pkin(1) + 0.2e1 * t1864) * t1768;
t1833 = (pkin(1) + 0.2e1 * t1863) * t1771;
t1832 = t1765 * t1766;
t1831 = t1766 * t1723;
t1830 = t1768 * t1769;
t1829 = t1769 * t1724;
t1828 = t1771 * t1772;
t1827 = t1772 * t1725;
t1826 = t1776 * t1786;
t1825 = t1779 * t1786;
t1824 = t1782 * t1786;
t1821 = t1765 * t1865;
t1820 = t1768 * t1864;
t1819 = t1771 * t1863;
t1687 = (-t1807 * t1705 + t1711 * t1714) * t1766 + (t1705 * t1714 + t1807 * t1711) * t1775;
t1729 = pkin(1) + t1865;
t1804 = pkin(2) * t1832 - t1729 * t1775;
t1818 = t1687 / t1804 * t1750;
t1688 = (-t1806 * t1706 + t1712 * t1715) * t1769 + (t1706 * t1715 + t1806 * t1712) * t1778;
t1731 = pkin(1) + t1864;
t1803 = pkin(2) * t1830 - t1731 * t1778;
t1817 = t1688 / t1803 * t1751;
t1689 = (-t1805 * t1707 + t1713 * t1716) * t1772 + (t1707 * t1716 + t1805 * t1713) * t1781;
t1733 = pkin(1) + t1863;
t1802 = pkin(2) * t1828 - t1733 * t1781;
t1816 = t1689 / t1802 * t1752;
t1815 = t1767 * t1832;
t1814 = t1770 * t1830;
t1813 = t1773 * t1828;
t1812 = t1776 * t1850;
t1811 = t1779 * t1849;
t1810 = t1782 * t1848;
t1801 = pkin(1) * t1815 + (t1815 * t1873 - t1826) * t1774;
t1800 = pkin(1) * t1814 + (t1814 * t1873 - t1825) * t1777;
t1799 = pkin(1) * t1813 + (t1813 * t1873 - t1824) * t1780;
t1798 = 0.1e1 / pkin(1);
t1794 = 0.2e1 * qJ(3,1);
t1791 = 0.2e1 * qJ(3,2);
t1788 = 0.2e1 * qJ(3,3);
t1758 = t1781 ^ 2;
t1756 = t1778 ^ 2;
t1754 = t1775 ^ 2;
t1704 = pkin(2) * t1771 * t1781 + t1733 * t1772;
t1703 = pkin(2) * t1768 * t1778 + t1731 * t1769;
t1702 = pkin(2) * t1765 * t1775 + t1729 * t1766;
t1698 = t1824 * t1828 + (t1757 - 0.1e1) * t1773 * pkin(2);
t1697 = t1825 * t1830 + (t1755 - 0.1e1) * t1770 * pkin(2);
t1696 = t1826 * t1832 + (t1753 - 0.1e1) * t1767 * pkin(2);
t1695 = -t1802 * t1773 + t1824;
t1694 = -t1803 * t1770 + t1825;
t1693 = -t1804 * t1767 + t1826;
t1 = [t1746 * t1812 + t1747 * t1811 + t1748 * t1810 - m(4) * g(1) + (-((t1745 * t1833 + t1748 * t1836) * t1758 + (t1745 * t1827 - t1799 * t1748) * t1781 - t1698 * t1748 - t1745 * t1819) * t1816 - ((t1744 * t1834 + t1747 * t1837) * t1756 + (t1744 * t1829 - t1800 * t1747) * t1778 - t1697 * t1747 - t1744 * t1820) * t1817 - ((t1743 * t1835 + t1746 * t1838) * t1754 + (t1743 * t1831 - t1801 * t1746) * t1775 - t1696 * t1746 - t1743 * t1821) * t1818 + ((-t1695 * t1748 - t1704 * t1745) * t1845 + (-t1694 * t1747 - t1703 * t1744) * t1846 + (-t1693 * t1746 - t1702 * t1743) * t1847) * t1866) * t1798; -t1743 * t1812 - t1744 * t1811 - t1745 * t1810 - m(4) * g(2) + (-((-t1745 * t1836 + t1748 * t1833) * t1758 + (t1799 * t1745 + t1748 * t1827) * t1781 + t1698 * t1745 - t1748 * t1819) * t1816 - ((-t1744 * t1837 + t1747 * t1834) * t1756 + (t1800 * t1744 + t1747 * t1829) * t1778 + t1697 * t1744 - t1747 * t1820) * t1817 - ((-t1743 * t1838 + t1746 * t1835) * t1754 + (t1801 * t1743 + t1746 * t1831) * t1775 + t1696 * t1743 - t1746 * t1821) * t1818 + ((t1695 * t1745 - t1704 * t1748) * t1845 + (t1694 * t1744 - t1703 * t1747) * t1846 + (t1693 * t1743 - t1702 * t1746) * t1847) * t1866) * t1798; -m(4) * g(3) - t1767 * t1850 - t1770 * t1849 - t1773 * t1848 + (((sin(-qJ(2,1) + t1862) + sin(qJ(2,1) + t1861)) * t1872 + (-cos(-0.2e1 * qJ(3,1) + t1855) - cos(t1794 + t1856) - 0.2e1 * t1782) * pkin(2) + (-cos(-qJ(3,1) + t1855) - cos(qJ(3,1) + t1856) - cos(t1862) - cos(t1861)) * pkin(1)) / ((-sin(t1794 + qJ(2,1)) + t1772) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t1739) * pkin(1)) * t1689 + ((sin(-qJ(2,2) + t1860) + sin(qJ(2,2) + t1859)) * t1872 + (-cos(-0.2e1 * qJ(3,2) + t1853) - cos(t1791 + t1854) - 0.2e1 * t1779) * pkin(2) + (-cos(-qJ(3,2) + t1853) - cos(qJ(3,2) + t1854) - cos(t1860) - cos(t1859)) * pkin(1)) / ((-sin(t1791 + qJ(2,2)) + t1769) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t1738) * pkin(1)) * t1688 + ((sin(-qJ(2,3) + t1858) + sin(qJ(2,3) + t1857)) * t1872 + (-cos(-0.2e1 * qJ(3,3) + t1851) - cos(t1788 + t1852) - 0.2e1 * t1776) * pkin(2) + (-cos(-qJ(3,3) + t1851) - cos(qJ(3,3) + t1852) - cos(t1858) - cos(t1857)) * pkin(1)) / ((-sin(t1788 + qJ(2,3)) + t1766) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t1737) * pkin(1)) * t1687) * t1798 / 0.2e1 + ((t1773 * t1786 + t1802 * t1782) * t1845 + (t1770 * t1786 + t1803 * t1779) * t1846 + (t1767 * t1786 + t1804 * t1776) * t1847) * t1798 * t1866;];
taugX  = t1;
