% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR12V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:47
% EndTime: 2020-08-06 18:21:48
% DurationCPUTime: 0.74s
% Computational Cost: add. (390->108), mult. (641->202), div. (63->7), fcn. (654->18), ass. (0->85)
t1838 = MDP(2) - MDP(4);
t1837 = -MDP(3) + MDP(5);
t1803 = sin(qJ(1,3));
t1809 = cos(qJ(1,3));
t1774 = -t1803 * g(1) + t1809 * g(2);
t1775 = t1809 * g(1) + t1803 * g(2);
t1802 = sin(qJ(3,3));
t1783 = t1802 * pkin(3) + qJ(2,3);
t1780 = 0.1e1 / t1783;
t1799 = legFrame(3,3);
t1786 = sin(t1799);
t1789 = cos(t1799);
t1832 = (-t1774 * t1786 - t1775 * t1789) * t1780;
t1805 = sin(qJ(1,2));
t1811 = cos(qJ(1,2));
t1776 = -t1805 * g(1) + t1811 * g(2);
t1777 = t1811 * g(1) + t1805 * g(2);
t1804 = sin(qJ(3,2));
t1784 = t1804 * pkin(3) + qJ(2,2);
t1781 = 0.1e1 / t1784;
t1800 = legFrame(2,3);
t1787 = sin(t1800);
t1790 = cos(t1800);
t1831 = (-t1776 * t1787 - t1777 * t1790) * t1781;
t1807 = sin(qJ(1,1));
t1813 = cos(qJ(1,1));
t1778 = g(1) * t1807 - g(2) * t1813;
t1779 = g(1) * t1813 + g(2) * t1807;
t1806 = sin(qJ(3,1));
t1785 = -t1806 * pkin(3) - qJ(2,1);
t1782 = 0.1e1 / t1785;
t1801 = legFrame(1,3);
t1788 = sin(t1801);
t1791 = cos(t1801);
t1830 = (t1778 * t1788 - t1779 * t1791) * t1782;
t1762 = -t1786 * t1803 + t1789 * t1809;
t1829 = t1762 * t1780;
t1763 = -t1787 * t1805 + t1790 * t1811;
t1828 = t1763 * t1781;
t1764 = -t1788 * t1807 + t1791 * t1813;
t1827 = t1764 * t1782;
t1765 = t1786 * t1809 + t1789 * t1803;
t1826 = t1765 * t1780;
t1766 = t1787 * t1811 + t1790 * t1805;
t1825 = t1766 * t1781;
t1767 = t1788 * t1813 + t1791 * t1807;
t1824 = t1767 * t1782;
t1823 = t1802 * t1832;
t1808 = cos(qJ(3,3));
t1822 = t1808 * t1832;
t1821 = t1804 * t1831;
t1810 = cos(qJ(3,2));
t1820 = t1810 * t1831;
t1819 = t1806 * t1830;
t1812 = cos(qJ(3,1));
t1818 = t1812 * t1830;
t1768 = t1786 * g(1) - t1789 * g(2);
t1771 = t1789 * g(1) + t1786 * g(2);
t1742 = t1768 * t1809 + t1771 * t1803;
t1741 = t1768 * t1803 - t1771 * t1809;
t1769 = t1787 * g(1) - t1790 * g(2);
t1772 = t1790 * g(1) + t1787 * g(2);
t1744 = t1769 * t1811 + t1772 * t1805;
t1743 = t1769 * t1805 - t1772 * t1811;
t1770 = t1788 * g(1) - t1791 * g(2);
t1773 = t1791 * g(1) + t1788 * g(2);
t1746 = t1770 * t1813 + t1773 * t1807;
t1745 = t1770 * t1807 - t1773 * t1813;
t1817 = t1774 * t1789 - t1775 * t1786;
t1816 = t1776 * t1790 - t1777 * t1787;
t1815 = t1778 * t1791 + t1779 * t1788;
t1798 = 0.1e1 / t1806;
t1797 = 0.1e1 / t1804;
t1796 = 0.1e1 / t1802;
t1795 = pkin(1) + pkin(5) + pkin(6);
t1761 = t1785 * t1813 + t1795 * t1807;
t1760 = -t1784 * t1811 + t1795 * t1805;
t1759 = -t1783 * t1809 + t1795 * t1803;
t1758 = -t1807 * t1785 + t1795 * t1813;
t1757 = t1805 * t1784 + t1795 * t1811;
t1756 = t1803 * t1783 + t1795 * t1809;
t1740 = -t1773 * (-t1807 * pkin(1) + t1813 * qJ(2,1)) + t1770 * (t1813 * pkin(1) + t1807 * qJ(2,1));
t1739 = -t1772 * (-t1805 * pkin(1) + t1811 * qJ(2,2)) + t1769 * (t1811 * pkin(1) + t1805 * qJ(2,2));
t1738 = -t1771 * (-t1803 * pkin(1) + t1809 * qJ(2,3)) + t1768 * (t1809 * pkin(1) + t1803 * qJ(2,3));
t1 = [(-(t1764 * t1740 - (t1758 * t1791 - t1788 * t1761) * t1746) * t1782 + (t1763 * t1739 - (t1757 * t1790 - t1787 * t1760) * t1744) * t1781 + (t1762 * t1738 - (t1756 * t1789 - t1786 * t1759) * t1742) * t1780) * MDP(6) + (t1762 * t1823 + t1763 * t1821 - t1764 * t1819) * MDP(12) + (t1762 * t1822 + t1763 * t1820 - t1764 * t1818) * MDP(13) - g(1) * MDP(14) + t1837 * (t1741 * t1829 + t1743 * t1828 - t1745 * t1827) + t1838 * (t1742 * t1829 + t1744 * t1828 - t1746 * t1827); (-(t1767 * t1740 - (t1788 * t1758 + t1761 * t1791) * t1746) * t1782 + (t1766 * t1739 - (t1787 * t1757 + t1760 * t1790) * t1744) * t1781 + (t1765 * t1738 - (t1786 * t1756 + t1759 * t1789) * t1742) * t1780) * MDP(6) + (t1765 * t1823 + t1766 * t1821 - t1767 * t1819) * MDP(12) + (t1765 * t1822 + t1766 * t1820 - t1767 * t1818) * MDP(13) - g(2) * MDP(14) + t1837 * (t1741 * t1826 + t1743 * t1825 - t1745 * t1824) + t1838 * (t1742 * t1826 + t1744 * t1825 - t1746 * t1824); (-t1742 * t1796 * t1808 - t1744 * t1797 * t1810 - t1746 * t1798 * t1812) * MDP(6) - g(3) * MDP(14) + ((-t1798 * (g(3) * t1806 - t1815 * t1812) - t1797 * (g(3) * t1804 + t1816 * t1810) - t1796 * (g(3) * t1802 + t1817 * t1808)) * MDP(12) + (-t1798 * (g(3) * t1812 + t1815 * t1806) - t1797 * (g(3) * t1810 - t1816 * t1804) - t1796 * (g(3) * t1808 - t1817 * t1802)) * MDP(13)) / pkin(3);];
taugX  = t1;
