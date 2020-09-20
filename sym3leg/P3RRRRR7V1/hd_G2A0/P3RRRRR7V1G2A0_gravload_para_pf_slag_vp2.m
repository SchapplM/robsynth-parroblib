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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
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

function taugX = P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:42:59
% EndTime: 2020-08-07 03:43:00
% DurationCPUTime: 1.27s
% Computational Cost: add. (855->236), mult. (1293->320), div. (60->14), fcn. (957->60), ass. (0->148)
t1814 = sin(qJ(3,1));
t1823 = cos(qJ(3,1));
t1901 = m(3) * pkin(1) + mrSges(2,1);
t1757 = mrSges(3,1) * t1823 - mrSges(3,2) * t1814 + t1901;
t1769 = t1814 * mrSges(3,1) + t1823 * mrSges(3,2) + mrSges(2,2);
t1815 = sin(qJ(2,1));
t1824 = cos(qJ(2,1));
t1910 = t1824 * t1757 - t1769 * t1815;
t1811 = sin(qJ(3,2));
t1820 = cos(qJ(3,2));
t1756 = mrSges(3,1) * t1820 - mrSges(3,2) * t1811 + t1901;
t1768 = t1811 * mrSges(3,1) + t1820 * mrSges(3,2) + mrSges(2,2);
t1812 = sin(qJ(2,2));
t1821 = cos(qJ(2,2));
t1909 = t1821 * t1756 - t1768 * t1812;
t1808 = sin(qJ(3,3));
t1817 = cos(qJ(3,3));
t1755 = mrSges(3,1) * t1817 - mrSges(3,2) * t1808 + t1901;
t1767 = t1808 * mrSges(3,1) + t1817 * mrSges(3,2) + mrSges(2,2);
t1809 = sin(qJ(2,3));
t1818 = cos(qJ(2,3));
t1908 = t1818 * t1755 - t1767 * t1809;
t1907 = 2 * pkin(2);
t1827 = (-pkin(5) - pkin(4));
t1906 = 2 * t1827;
t1904 = t1817 * pkin(2);
t1903 = t1820 * pkin(2);
t1902 = t1823 * pkin(2);
t1900 = -qJ(3,1) + qJ(1,1);
t1899 = qJ(3,1) + qJ(1,1);
t1898 = -qJ(3,2) + qJ(1,2);
t1897 = qJ(3,2) + qJ(1,2);
t1896 = -qJ(3,3) + qJ(1,3);
t1895 = qJ(3,3) + qJ(1,3);
t1894 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1893 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1892 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1891 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1890 = qJ(1,3) + 0.2e1 * qJ(2,3);
t1889 = -0.2e1 * qJ(2,3) + qJ(1,3);
t1805 = legFrame(3,2);
t1787 = sin(t1805);
t1790 = cos(t1805);
t1758 = t1787 * g(1) + t1790 * g(2);
t1802 = qJ(2,3) + qJ(3,3);
t1781 = sin(t1802);
t1784 = cos(t1802);
t1793 = 0.1e1 / t1808;
t1761 = t1790 * g(1) - t1787 * g(2);
t1810 = sin(qJ(1,3));
t1819 = cos(qJ(1,3));
t1848 = g(3) * t1819 + t1761 * t1810;
t1888 = ((-mrSges(3,1) * t1758 + t1848 * mrSges(3,2)) * t1784 + t1781 * (t1848 * mrSges(3,1) + mrSges(3,2) * t1758)) * t1793;
t1806 = legFrame(2,2);
t1788 = sin(t1806);
t1791 = cos(t1806);
t1759 = t1788 * g(1) + t1791 * g(2);
t1803 = qJ(2,2) + qJ(3,2);
t1782 = sin(t1803);
t1785 = cos(t1803);
t1794 = 0.1e1 / t1811;
t1762 = t1791 * g(1) - t1788 * g(2);
t1813 = sin(qJ(1,2));
t1822 = cos(qJ(1,2));
t1847 = g(3) * t1822 + t1762 * t1813;
t1887 = ((-mrSges(3,1) * t1759 + t1847 * mrSges(3,2)) * t1785 + t1782 * (t1847 * mrSges(3,1) + mrSges(3,2) * t1759)) * t1794;
t1807 = legFrame(1,2);
t1789 = sin(t1807);
t1792 = cos(t1807);
t1760 = t1789 * g(1) + t1792 * g(2);
t1804 = qJ(2,1) + qJ(3,1);
t1783 = sin(t1804);
t1786 = cos(t1804);
t1795 = 0.1e1 / t1814;
t1763 = t1792 * g(1) - t1789 * g(2);
t1816 = sin(qJ(1,1));
t1825 = cos(qJ(1,1));
t1846 = g(3) * t1825 + t1763 * t1816;
t1886 = ((-mrSges(3,1) * t1760 + t1846 * mrSges(3,2)) * t1786 + t1783 * (t1846 * mrSges(3,1) + mrSges(3,2) * t1760)) * t1795;
t1774 = m(3) * pkin(4) - mrSges(1,2) + mrSges(2,3) + mrSges(3,3);
t1773 = t1774 * g(3);
t1826 = mrSges(1,1) * g(3);
t1885 = 0.1e1 / (t1818 * pkin(1) + pkin(2) * t1784) * (-t1773 * t1819 + t1810 * (t1908 * g(3) + t1826) + ((-mrSges(1,1) - t1908) * t1819 - t1810 * t1774) * t1761);
t1884 = 0.1e1 / (t1821 * pkin(1) + pkin(2) * t1785) * (-t1773 * t1822 + t1813 * (t1909 * g(3) + t1826) + ((-mrSges(1,1) - t1909) * t1822 - t1813 * t1774) * t1762);
t1883 = 0.1e1 / (t1824 * pkin(1) + pkin(2) * t1786) * (-t1773 * t1825 + t1816 * (t1910 * g(3) + t1826) + ((-mrSges(1,1) - t1910) * t1825 - t1816 * t1774) * t1763);
t1796 = t1817 ^ 2;
t1770 = pkin(1) * t1817 + t1796 * t1907 - pkin(2);
t1879 = t1770 * t1810;
t1798 = t1820 ^ 2;
t1771 = pkin(1) * t1820 + t1798 * t1907 - pkin(2);
t1878 = t1771 * t1813;
t1800 = t1823 ^ 2;
t1772 = pkin(1) * t1823 + t1800 * t1907 - pkin(2);
t1877 = t1772 * t1816;
t1876 = (pkin(1) + 0.2e1 * t1904) * t1808;
t1875 = (pkin(1) + 0.2e1 * t1903) * t1811;
t1874 = (pkin(1) + 0.2e1 * t1902) * t1814;
t1873 = t1808 * t1809;
t1872 = t1809 * t1770;
t1871 = t1811 * t1812;
t1870 = t1812 * t1771;
t1869 = t1814 * t1815;
t1868 = t1815 * t1772;
t1866 = t1819 * t1827;
t1864 = t1822 * t1827;
t1862 = t1825 * t1827;
t1860 = t1808 * t1904;
t1859 = t1811 * t1903;
t1858 = t1814 * t1902;
t1731 = (t1848 * t1755 + t1758 * t1767) * t1809 - t1818 * (t1758 * t1755 - t1848 * t1767);
t1776 = pkin(1) + t1904;
t1845 = pkin(2) * t1873 - t1776 * t1818;
t1857 = t1731 / t1845 * t1793;
t1732 = (t1847 * t1756 + t1759 * t1768) * t1812 - t1821 * (t1759 * t1756 - t1847 * t1768);
t1778 = pkin(1) + t1903;
t1844 = pkin(2) * t1871 - t1778 * t1821;
t1856 = t1732 / t1844 * t1794;
t1733 = (t1846 * t1757 + t1760 * t1769) * t1815 - t1824 * (t1760 * t1757 - t1846 * t1769);
t1780 = pkin(1) + t1902;
t1843 = pkin(2) * t1869 - t1780 * t1824;
t1855 = t1733 / t1843 * t1795;
t1854 = t1810 * t1873;
t1853 = t1813 * t1871;
t1852 = t1816 * t1869;
t1851 = t1819 * t1885;
t1850 = t1822 * t1884;
t1849 = t1825 * t1883;
t1842 = pkin(1) * t1854 + (t1854 * t1907 - t1866) * t1817;
t1841 = pkin(1) * t1853 + (t1853 * t1907 - t1864) * t1820;
t1840 = pkin(1) * t1852 + (t1852 * t1907 - t1862) * t1823;
t1839 = 0.1e1 / pkin(1);
t1838 = 1 / pkin(2);
t1835 = 0.2e1 * qJ(3,1);
t1832 = 0.2e1 * qJ(3,2);
t1829 = 0.2e1 * qJ(3,3);
t1801 = t1824 ^ 2;
t1799 = t1821 ^ 2;
t1797 = t1818 ^ 2;
t1751 = t1814 * pkin(2) * t1824 + t1815 * t1780;
t1750 = t1811 * pkin(2) * t1821 + t1812 * t1778;
t1749 = t1808 * pkin(2) * t1818 + t1809 * t1776;
t1745 = t1862 * t1869 + (t1800 - 0.1e1) * t1816 * pkin(2);
t1744 = t1864 * t1871 + (t1798 - 0.1e1) * t1813 * pkin(2);
t1743 = t1866 * t1873 + (t1796 - 0.1e1) * t1810 * pkin(2);
t1742 = -t1843 * t1816 + t1862;
t1741 = -t1844 * t1813 + t1864;
t1740 = -t1845 * t1810 + t1866;
t1 = [t1790 * t1851 + t1791 * t1850 + t1792 * t1849 - g(1) * m(4) + (-((t1789 * t1874 + t1792 * t1877) * t1801 + (t1789 * t1868 - t1840 * t1792) * t1824 - t1745 * t1792 - t1789 * t1858) * t1855 - ((t1788 * t1875 + t1791 * t1878) * t1799 + (t1788 * t1870 - t1841 * t1791) * t1821 - t1744 * t1791 - t1788 * t1859) * t1856 - ((t1787 * t1876 + t1790 * t1879) * t1797 + (t1787 * t1872 - t1842 * t1790) * t1818 - t1743 * t1790 - t1787 * t1860) * t1857 + ((-t1742 * t1792 - t1751 * t1789) * t1886 + (-t1741 * t1791 - t1750 * t1788) * t1887 + (-t1740 * t1790 - t1749 * t1787) * t1888) * t1838) * t1839; -t1787 * t1851 - t1788 * t1850 - t1789 * t1849 - g(2) * m(4) + (-((-t1789 * t1877 + t1792 * t1874) * t1801 + (t1840 * t1789 + t1792 * t1868) * t1824 + t1745 * t1789 - t1792 * t1858) * t1855 - ((-t1788 * t1878 + t1791 * t1875) * t1799 + (t1841 * t1788 + t1791 * t1870) * t1821 + t1744 * t1788 - t1791 * t1859) * t1856 - ((-t1787 * t1879 + t1790 * t1876) * t1797 + (t1842 * t1787 + t1790 * t1872) * t1818 + t1743 * t1787 - t1790 * t1860) * t1857 + ((t1742 * t1789 - t1751 * t1792) * t1886 + (t1741 * t1788 - t1750 * t1791) * t1887 + (t1740 * t1787 - t1749 * t1790) * t1888) * t1838) * t1839; -g(3) * m(4) - t1810 * t1885 - t1813 * t1884 - t1816 * t1883 + (((sin(-qJ(2,1) + t1900) + sin(qJ(2,1) + t1899)) * t1906 + (-cos(-0.2e1 * qJ(3,1) + t1893) - cos(t1835 + t1894) - 0.2e1 * t1825) * pkin(2) + (-cos(-qJ(3,1) + t1893) - cos(qJ(3,1) + t1894) - cos(t1900) - cos(t1899)) * pkin(1)) / ((-sin(t1835 + qJ(2,1)) + t1815) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t1783) * pkin(1)) * t1733 + ((sin(-qJ(2,2) + t1898) + sin(qJ(2,2) + t1897)) * t1906 + (-cos(-0.2e1 * qJ(3,2) + t1891) - cos(t1832 + t1892) - 0.2e1 * t1822) * pkin(2) + (-cos(-qJ(3,2) + t1891) - cos(qJ(3,2) + t1892) - cos(t1898) - cos(t1897)) * pkin(1)) / ((-sin(t1832 + qJ(2,2)) + t1812) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t1782) * pkin(1)) * t1732 + ((sin(-qJ(2,3) + t1896) + sin(qJ(2,3) + t1895)) * t1906 + (-cos(-0.2e1 * qJ(3,3) + t1889) - cos(t1829 + t1890) - 0.2e1 * t1819) * pkin(2) + (-cos(-qJ(3,3) + t1889) - cos(qJ(3,3) + t1890) - cos(t1896) - cos(t1895)) * pkin(1)) / ((-sin(t1829 + qJ(2,3)) + t1809) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t1781) * pkin(1)) * t1731) * t1839 / 0.2e1 + ((t1816 * t1827 + t1843 * t1825) * t1886 + (t1813 * t1827 + t1844 * t1822) * t1887 + (t1810 * t1827 + t1845 * t1819) * t1888) * t1838 * t1839;];
taugX  = t1;
