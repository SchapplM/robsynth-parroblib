% Calculate Gravitation load for parallel robot
% P3RRRRR7V1G3A0
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
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:52:10
% EndTime: 2020-08-07 08:52:12
% DurationCPUTime: 1.60s
% Computational Cost: add. (864->236), mult. (1656->330), div. (60->14), fcn. (969->60), ass. (0->152)
t1821 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - m(3) * (rSges(3,3) + pkin(4));
t1779 = sin(qJ(1,3));
t1798 = -pkin(5) - pkin(4);
t1742 = t1779 * t1798;
t1788 = cos(qJ(1,3));
t1786 = cos(qJ(3,3));
t1874 = pkin(2) * t1786;
t1737 = pkin(1) + t1874;
t1787 = cos(qJ(2,3));
t1777 = sin(qJ(3,3));
t1778 = sin(qJ(2,3));
t1841 = t1777 * t1778;
t1816 = pkin(2) * t1841 - t1737 * t1787;
t1888 = t1816 * t1788 + t1742;
t1782 = sin(qJ(1,2));
t1743 = t1782 * t1798;
t1791 = cos(qJ(1,2));
t1789 = cos(qJ(3,2));
t1873 = pkin(2) * t1789;
t1739 = pkin(1) + t1873;
t1790 = cos(qJ(2,2));
t1780 = sin(qJ(3,2));
t1781 = sin(qJ(2,2));
t1839 = t1780 * t1781;
t1815 = pkin(2) * t1839 - t1739 * t1790;
t1887 = t1815 * t1791 + t1743;
t1785 = sin(qJ(1,1));
t1744 = t1785 * t1798;
t1794 = cos(qJ(1,1));
t1792 = cos(qJ(3,1));
t1872 = pkin(2) * t1792;
t1741 = pkin(1) + t1872;
t1793 = cos(qJ(2,1));
t1783 = sin(qJ(3,1));
t1784 = sin(qJ(2,1));
t1837 = t1783 * t1784;
t1814 = pkin(2) * t1837 - t1741 * t1793;
t1886 = t1814 * t1794 + t1744;
t1835 = -m(2) * rSges(2,1) - pkin(1) * m(3);
t1878 = m(3) * rSges(3,2);
t1879 = m(3) * rSges(3,1);
t1716 = t1783 * t1878 - t1792 * t1879 + t1835;
t1796 = m(2) * rSges(2,2);
t1722 = t1796 + (rSges(3,1) * t1783 + rSges(3,2) * t1792) * m(3);
t1885 = -t1716 * t1793 - t1722 * t1784;
t1715 = t1780 * t1878 - t1789 * t1879 + t1835;
t1721 = t1796 + (rSges(3,1) * t1780 + rSges(3,2) * t1789) * m(3);
t1884 = -t1715 * t1790 - t1721 * t1781;
t1714 = t1777 * t1878 - t1786 * t1879 + t1835;
t1720 = t1796 + (rSges(3,1) * t1777 + rSges(3,2) * t1786) * m(3);
t1883 = -t1714 * t1787 - t1720 * t1778;
t1882 = 0.2e1 * pkin(2);
t1881 = 0.2e1 * t1798;
t1880 = m(1) * rSges(1,1);
t1875 = m(3) / pkin(2);
t1871 = -qJ(3,1) + qJ(1,1);
t1870 = qJ(3,1) + qJ(1,1);
t1869 = -qJ(3,2) + qJ(1,2);
t1868 = qJ(3,2) + qJ(1,2);
t1867 = -qJ(3,3) + qJ(1,3);
t1866 = qJ(3,3) + qJ(1,3);
t1865 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1864 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1863 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1862 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1861 = qJ(1,3) + 0.2e1 * qJ(2,3);
t1860 = qJ(1,3) - 0.2e1 * qJ(2,3);
t1774 = legFrame(3,2);
t1755 = sin(t1774);
t1758 = cos(t1774);
t1723 = g(1) * t1755 + g(2) * t1758;
t1771 = qJ(2,3) + qJ(3,3);
t1749 = sin(t1771);
t1752 = cos(t1771);
t1762 = 0.1e1 / t1777;
t1726 = g(1) * t1758 - g(2) * t1755;
t1819 = g(3) * t1779 - t1726 * t1788;
t1859 = ((-rSges(3,1) * t1723 - t1819 * rSges(3,2)) * t1752 + t1749 * (-t1819 * rSges(3,1) + rSges(3,2) * t1723)) * t1762;
t1775 = legFrame(2,2);
t1756 = sin(t1775);
t1759 = cos(t1775);
t1724 = g(1) * t1756 + g(2) * t1759;
t1772 = qJ(2,2) + qJ(3,2);
t1750 = sin(t1772);
t1753 = cos(t1772);
t1763 = 0.1e1 / t1780;
t1727 = g(1) * t1759 - g(2) * t1756;
t1818 = g(3) * t1782 - t1727 * t1791;
t1858 = ((-rSges(3,1) * t1724 - t1818 * rSges(3,2)) * t1753 + t1750 * (-t1818 * rSges(3,1) + rSges(3,2) * t1724)) * t1763;
t1776 = legFrame(1,2);
t1757 = sin(t1776);
t1760 = cos(t1776);
t1725 = g(1) * t1757 + g(2) * t1760;
t1773 = qJ(2,1) + qJ(3,1);
t1751 = sin(t1773);
t1754 = cos(t1773);
t1764 = 0.1e1 / t1783;
t1728 = g(1) * t1760 - g(2) * t1757;
t1817 = g(3) * t1785 - t1728 * t1794;
t1857 = ((-rSges(3,1) * t1725 - t1817 * rSges(3,2)) * t1754 + t1751 * (-t1817 * rSges(3,1) + rSges(3,2) * t1725)) * t1764;
t1761 = g(3) * t1880;
t1820 = t1821 * g(3);
t1850 = 0.1e1 / (pkin(1) * t1787 + pkin(2) * t1752) * ((g(3) * t1883 + t1761) * t1788 - t1820 * t1779 + (t1821 * t1788 + (t1880 + t1883) * t1779) * t1726);
t1849 = 0.1e1 / (pkin(1) * t1790 + pkin(2) * t1753) * ((g(3) * t1884 + t1761) * t1791 - t1820 * t1782 + (t1821 * t1791 + (t1880 + t1884) * t1782) * t1727);
t1848 = 0.1e1 / (pkin(1) * t1793 + pkin(2) * t1754) * ((g(3) * t1885 + t1761) * t1794 - t1820 * t1785 + (t1821 * t1794 + (t1880 + t1885) * t1785) * t1728);
t1765 = t1786 ^ 2;
t1732 = pkin(1) * t1786 + t1765 * t1882 - pkin(2);
t1847 = t1732 * t1788;
t1767 = t1789 ^ 2;
t1733 = pkin(1) * t1789 + t1767 * t1882 - pkin(2);
t1846 = t1733 * t1791;
t1769 = t1792 ^ 2;
t1734 = pkin(1) * t1792 + t1769 * t1882 - pkin(2);
t1845 = t1734 * t1794;
t1844 = (pkin(1) + 0.2e1 * t1874) * t1777;
t1843 = (pkin(1) + 0.2e1 * t1873) * t1780;
t1842 = (pkin(1) + 0.2e1 * t1872) * t1783;
t1840 = t1778 * t1732;
t1838 = t1781 * t1733;
t1836 = t1784 * t1734;
t1833 = t1777 * t1874;
t1832 = t1780 * t1873;
t1831 = t1783 * t1872;
t1696 = (t1819 * t1714 + t1720 * t1723) * t1778 + (t1714 * t1723 - t1819 * t1720) * t1787;
t1830 = t1696 / t1816 * t1762;
t1697 = (t1818 * t1715 + t1721 * t1724) * t1781 + (t1715 * t1724 - t1818 * t1721) * t1790;
t1829 = t1697 / t1815 * t1763;
t1698 = (t1817 * t1716 + t1722 * t1725) * t1784 + (t1716 * t1725 - t1817 * t1722) * t1793;
t1828 = t1698 / t1814 * t1764;
t1827 = t1788 * t1841;
t1826 = t1779 * t1850;
t1825 = t1791 * t1839;
t1824 = t1782 * t1849;
t1823 = t1794 * t1837;
t1822 = t1785 * t1848;
t1813 = pkin(1) * t1827 + (t1827 * t1882 + t1742) * t1786;
t1812 = pkin(1) * t1825 + (t1825 * t1882 + t1743) * t1789;
t1811 = pkin(1) * t1823 + (t1823 * t1882 + t1744) * t1792;
t1810 = 0.1e1 / pkin(1);
t1806 = 0.2e1 * qJ(3,1);
t1803 = 0.2e1 * qJ(3,2);
t1800 = 0.2e1 * qJ(3,3);
t1770 = t1793 ^ 2;
t1768 = t1790 ^ 2;
t1766 = t1787 ^ 2;
t1710 = pkin(2) * t1783 * t1793 + t1741 * t1784;
t1709 = pkin(2) * t1780 * t1790 + t1739 * t1781;
t1708 = pkin(2) * t1777 * t1787 + t1737 * t1778;
t1704 = -t1744 * t1837 + (t1769 - 0.1e1) * t1794 * pkin(2);
t1703 = -t1743 * t1839 + (t1767 - 0.1e1) * t1791 * pkin(2);
t1702 = -t1742 * t1841 + (t1765 - 0.1e1) * t1788 * pkin(2);
t1 = [-t1758 * t1826 - t1759 * t1824 - t1760 * t1822 - m(4) * g(1) + (-((t1757 * t1842 + t1760 * t1845) * t1770 + (t1757 * t1836 - t1811 * t1760) * t1793 - t1704 * t1760 - t1757 * t1831) * t1828 - ((t1756 * t1843 + t1759 * t1846) * t1768 + (t1756 * t1838 - t1812 * t1759) * t1790 - t1703 * t1759 - t1756 * t1832) * t1829 - ((t1755 * t1844 + t1758 * t1847) * t1766 + (t1755 * t1840 - t1813 * t1758) * t1787 - t1702 * t1758 - t1755 * t1833) * t1830 + ((-t1757 * t1710 + t1886 * t1760) * t1857 + (-t1756 * t1709 + t1887 * t1759) * t1858 + (-t1755 * t1708 + t1888 * t1758) * t1859) * t1875) * t1810; t1755 * t1826 + t1756 * t1824 + t1757 * t1822 - m(4) * g(2) + (-((-t1757 * t1845 + t1760 * t1842) * t1770 + (t1811 * t1757 + t1760 * t1836) * t1793 + t1704 * t1757 - t1760 * t1831) * t1828 - ((-t1756 * t1846 + t1759 * t1843) * t1768 + (t1812 * t1756 + t1759 * t1838) * t1790 + t1703 * t1756 - t1759 * t1832) * t1829 - ((-t1755 * t1847 + t1758 * t1844) * t1766 + (t1813 * t1755 + t1758 * t1840) * t1787 + t1702 * t1755 - t1758 * t1833) * t1830 + ((-t1760 * t1710 - t1886 * t1757) * t1857 + (-t1759 * t1709 - t1887 * t1756) * t1858 + (-t1758 * t1708 - t1888 * t1755) * t1859) * t1875) * t1810; -m(4) * g(3) - t1788 * t1850 - t1791 * t1849 - t1794 * t1848 + (((cos(-qJ(2,1) + t1871) + cos(qJ(2,1) + t1870)) * t1881 + (sin(-0.2e1 * qJ(3,1) + t1864) + sin(t1806 + t1865) + 0.2e1 * t1785) * pkin(2) + (sin(-qJ(3,1) + t1864) + sin(qJ(3,1) + t1865) + sin(t1871) + sin(t1870)) * pkin(1)) / ((-sin(t1806 + qJ(2,1)) + t1784) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t1751) * pkin(1)) * t1698 + ((cos(-qJ(2,2) + t1869) + cos(qJ(2,2) + t1868)) * t1881 + (sin(-0.2e1 * qJ(3,2) + t1862) + sin(t1803 + t1863) + 0.2e1 * t1782) * pkin(2) + (sin(-qJ(3,2) + t1862) + sin(qJ(3,2) + t1863) + sin(t1869) + sin(t1868)) * pkin(1)) / ((-sin(t1803 + qJ(2,2)) + t1781) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t1750) * pkin(1)) * t1697 + ((cos(-qJ(2,3) + t1867) + cos(qJ(2,3) + t1866)) * t1881 + (sin(-0.2e1 * qJ(3,3) + t1860) + sin(t1800 + t1861) + 0.2e1 * t1779) * pkin(2) + (sin(-qJ(3,3) + t1860) + sin(qJ(3,3) + t1861) + sin(t1867) + sin(t1866)) * pkin(1)) / ((-sin(t1800 + qJ(2,3)) + t1778) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t1749) * pkin(1)) * t1696) * t1810 / 0.2e1 + ((-t1785 * t1814 + t1794 * t1798) * t1857 + (-t1782 * t1815 + t1791 * t1798) * t1858 + (-t1779 * t1816 + t1788 * t1798) * t1859) * t1810 * t1875;];
taugX  = t1;
