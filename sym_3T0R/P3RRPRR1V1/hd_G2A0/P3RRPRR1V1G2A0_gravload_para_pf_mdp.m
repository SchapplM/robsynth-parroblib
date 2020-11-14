% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR1V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR1V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:34:18
% EndTime: 2020-08-06 19:34:19
% DurationCPUTime: 0.95s
% Computational Cost: add. (453->117), mult. (759->223), div. (129->7), fcn. (789->18), ass. (0->104)
t1957 = MDP(3) - MDP(11);
t1884 = sin(qJ(1,3));
t1890 = cos(qJ(1,3));
t1880 = legFrame(3,2);
t1864 = sin(t1880);
t1867 = cos(t1880);
t1909 = t1867 * g(1) - t1864 * g(2);
t1846 = g(3) * t1884 - t1890 * t1909;
t1877 = pkin(3) + qJ(3,3);
t1870 = 0.1e1 / t1877;
t1956 = t1846 * t1870;
t1938 = t1870 * t1890;
t1924 = t1846 * t1938;
t1886 = sin(qJ(1,2));
t1892 = cos(qJ(1,2));
t1881 = legFrame(2,2);
t1865 = sin(t1881);
t1868 = cos(t1881);
t1908 = t1868 * g(1) - t1865 * g(2);
t1847 = g(3) * t1886 - t1892 * t1908;
t1878 = pkin(3) + qJ(3,2);
t1871 = 0.1e1 / t1878;
t1955 = t1847 * t1871;
t1936 = t1871 * t1892;
t1923 = t1847 * t1936;
t1888 = sin(qJ(1,1));
t1894 = cos(qJ(1,1));
t1882 = legFrame(1,2);
t1866 = sin(t1882);
t1869 = cos(t1882);
t1907 = t1869 * g(1) - t1866 * g(2);
t1848 = g(3) * t1888 - t1894 * t1907;
t1879 = pkin(3) + qJ(3,1);
t1872 = 0.1e1 / t1879;
t1954 = t1848 * t1872;
t1934 = t1872 * t1894;
t1922 = t1848 * t1934;
t1953 = MDP(12) * pkin(1) + MDP(9);
t1889 = cos(qJ(2,3));
t1873 = 0.1e1 / t1889;
t1945 = t1864 * t1873;
t1891 = cos(qJ(2,2));
t1874 = 0.1e1 / t1891;
t1944 = t1865 * t1874;
t1893 = cos(qJ(2,1));
t1875 = 0.1e1 / t1893;
t1943 = t1866 * t1875;
t1942 = t1867 * t1873;
t1941 = t1868 * t1874;
t1940 = t1869 * t1875;
t1939 = t1870 * t1873;
t1937 = t1871 * t1874;
t1935 = t1872 * t1875;
t1883 = sin(qJ(2,3));
t1895 = pkin(1) + pkin(2);
t1933 = t1883 * t1895;
t1932 = t1884 * t1889;
t1885 = sin(qJ(2,2));
t1931 = t1885 * t1895;
t1930 = t1886 * t1891;
t1887 = sin(qJ(2,1));
t1929 = t1887 * t1895;
t1928 = t1888 * t1893;
t1927 = t1889 * t1895;
t1926 = t1891 * t1895;
t1925 = t1893 * t1895;
t1852 = -t1864 * t1932 + t1883 * t1867;
t1921 = t1852 * t1939;
t1853 = -t1865 * t1930 + t1885 * t1868;
t1920 = t1853 * t1937;
t1854 = -t1866 * t1928 + t1887 * t1869;
t1919 = t1854 * t1935;
t1855 = t1883 * t1864 + t1867 * t1932;
t1918 = t1855 * t1939;
t1856 = t1885 * t1865 + t1868 * t1930;
t1917 = t1856 * t1937;
t1857 = t1887 * t1866 + t1869 * t1928;
t1916 = t1857 * t1935;
t1837 = g(3) * (pkin(1) * t1932 - t1890 * qJ(3,3)) - t1909 * (t1890 * t1889 * pkin(1) + t1884 * qJ(3,3));
t1915 = t1837 * t1939;
t1838 = g(3) * (pkin(1) * t1930 - t1892 * qJ(3,2)) - t1908 * (t1892 * t1891 * pkin(1) + t1886 * qJ(3,2));
t1914 = t1838 * t1937;
t1839 = g(3) * (pkin(1) * t1928 - t1894 * qJ(3,1)) - t1907 * (t1894 * t1893 * pkin(1) + t1888 * qJ(3,1));
t1913 = t1839 * t1935;
t1912 = t1846 * t1883 * t1939;
t1911 = t1847 * t1885 * t1937;
t1910 = t1848 * t1887 * t1935;
t1906 = -t1877 * t1890 + t1884 * t1927;
t1905 = -t1878 * t1892 + t1886 * t1926;
t1904 = -t1879 * t1894 + t1888 * t1925;
t1903 = g(3) * t1890 + t1884 * t1909;
t1902 = g(3) * t1892 + t1886 * t1908;
t1901 = g(3) * t1894 + t1888 * t1907;
t1876 = 0.1e1 / t1895;
t1863 = t1866 * g(1) + t1869 * g(2);
t1862 = t1865 * g(1) + t1868 * g(2);
t1861 = t1864 * g(1) + t1867 * g(2);
t1836 = t1863 * t1887 + t1893 * t1901;
t1835 = -t1863 * t1893 + t1887 * t1901;
t1834 = t1862 * t1885 + t1891 * t1902;
t1833 = -t1862 * t1891 + t1885 * t1902;
t1832 = t1861 * t1883 + t1889 * t1903;
t1831 = -t1861 * t1889 + t1883 * t1903;
t1 = [(t1846 * t1918 + t1847 * t1917 + t1848 * t1916) * MDP(2) + (t1855 * t1956 + t1856 * t1955 + t1857 * t1954) * MDP(9) + (-t1855 * t1912 - t1856 * t1911 - t1857 * t1910) * MDP(10) + (t1857 * t1913 - (t1866 * t1929 + t1869 * t1904) * t1954 + t1856 * t1914 - (t1865 * t1931 + t1868 * t1905) * t1955 + t1855 * t1915 - (t1864 * t1933 + t1867 * t1906) * t1956) * MDP(12) - g(1) * MDP(13) + ((t1832 * t1945 + t1834 * t1944 + t1836 * t1943) * MDP(10) + t1953 * (t1831 * t1945 + t1833 * t1944 + t1835 * t1943)) * t1876 + t1957 * (t1901 * t1916 + t1902 * t1917 + t1903 * t1918); (t1846 * t1921 + t1847 * t1920 + t1848 * t1919) * MDP(2) + (t1852 * t1956 + t1853 * t1955 + t1854 * t1954) * MDP(9) + (-t1852 * t1912 - t1853 * t1911 - t1854 * t1910) * MDP(10) + (t1854 * t1913 - (-t1866 * t1904 + t1869 * t1929) * t1954 + t1853 * t1914 - (-t1865 * t1905 + t1868 * t1931) * t1955 + t1852 * t1915 - (-t1864 * t1906 + t1867 * t1933) * t1956) * MDP(12) - g(2) * MDP(13) + ((t1832 * t1942 + t1834 * t1941 + t1836 * t1940) * MDP(10) + t1953 * (t1831 * t1942 + t1833 * t1941 + t1835 * t1940)) * t1876 + t1957 * (t1901 * t1919 + t1902 * t1920 + t1903 * t1921); (t1922 + t1923 + t1924) * MDP(2) + (t1889 * t1924 + t1891 * t1923 + t1893 * t1922) * MDP(9) + (-t1883 * t1924 - t1885 * t1923 - t1887 * t1922) * MDP(10) + ((t1894 * t1839 - (t1888 * t1879 + t1894 * t1925) * t1848) * t1872 + (t1892 * t1838 - (t1886 * t1878 + t1892 * t1926) * t1847) * t1871 + (t1890 * t1837 - (t1884 * t1877 + t1890 * t1927) * t1846) * t1870) * MDP(12) - g(3) * MDP(13) + t1957 * (t1901 * t1934 + t1902 * t1936 + t1903 * t1938);];
taugX  = t1;
