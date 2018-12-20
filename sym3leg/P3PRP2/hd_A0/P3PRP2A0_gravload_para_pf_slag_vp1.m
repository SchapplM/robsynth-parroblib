% Calculate Gravitation load for parallel robot
% P3PRP2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (including platform)
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
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3PRP2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:38:55
% EndTime: 2018-12-20 17:38:56
% DurationCPUTime: 0.53s
% Computational Cost: add. (690->153), mult. (1242->254), div. (36->3), fcn. (644->14), ass. (0->124)
t822 = (pkin(2) ^ 2);
t856 = -t822 - 1;
t802 = legFrame(1,3);
t787 = sin(t802);
t855 = qJ(3,1) * t787;
t790 = cos(t802);
t854 = qJ(3,1) * t790;
t801 = legFrame(2,3);
t786 = sin(t801);
t853 = qJ(3,2) * t786;
t789 = cos(t801);
t852 = qJ(3,2) * t789;
t800 = legFrame(3,3);
t785 = sin(t800);
t851 = qJ(3,3) * t785;
t788 = cos(t800);
t850 = qJ(3,3) * t788;
t770 = -t785 * g(1) + t788 * g(2);
t773 = t788 * g(1) + t785 * g(2);
t797 = rSges(3,3) + qJ(3,3);
t803 = sin(qJ(2,3));
t806 = cos(qJ(2,3));
t809 = pkin(2) + rSges(3,1);
t731 = ((-t770 * t797 - t809 * t773) * m(3) - m(2) * (rSges(2,1) * t773 - rSges(2,2) * t770)) * t806 + ((t770 * t809 - t797 * t773) * m(3) + m(2) * (rSges(2,1) * t770 + rSges(2,2) * t773)) * t803;
t813 = qJ(3,3) ^ 2;
t782 = -t813 - t856;
t794 = t806 ^ 2;
t837 = t806 * qJ(3,3);
t831 = 0.2e1 * t837;
t755 = 0.1e1 / (t803 * pkin(2) * t831 + t782 * t794 - t813 + t856);
t849 = t731 * t755;
t771 = -t786 * g(1) + t789 * g(2);
t774 = t789 * g(1) + t786 * g(2);
t798 = rSges(3,3) + qJ(3,2);
t804 = sin(qJ(2,2));
t807 = cos(qJ(2,2));
t732 = ((-t771 * t798 - t809 * t774) * m(3) - m(2) * (rSges(2,1) * t774 - rSges(2,2) * t771)) * t807 + ((t771 * t809 - t798 * t774) * m(3) + m(2) * (rSges(2,1) * t771 + rSges(2,2) * t774)) * t804;
t814 = qJ(3,2) ^ 2;
t783 = -t814 - t856;
t795 = t807 ^ 2;
t836 = t807 * qJ(3,2);
t830 = 0.2e1 * t836;
t756 = 0.1e1 / (t804 * pkin(2) * t830 + t783 * t795 - t814 + t856);
t848 = t732 * t756;
t772 = -t787 * g(1) + t790 * g(2);
t775 = t790 * g(1) + t787 * g(2);
t799 = rSges(3,3) + qJ(3,1);
t805 = sin(qJ(2,1));
t808 = cos(qJ(2,1));
t733 = ((-t772 * t799 - t809 * t775) * m(3) - m(2) * (rSges(2,1) * t775 - rSges(2,2) * t772)) * t808 + ((t772 * t809 - t799 * t775) * m(3) + m(2) * (rSges(2,1) * t772 + rSges(2,2) * t775)) * t805;
t815 = qJ(3,1) ^ 2;
t784 = -t815 - t856;
t796 = t808 ^ 2;
t835 = t808 * qJ(3,1);
t829 = 0.2e1 * t835;
t757 = 0.1e1 / (t805 * pkin(2) * t829 + t784 * t796 - t815 + t856);
t847 = t733 * t757;
t746 = -t770 * t803 + t773 * t806;
t846 = t746 * t755;
t747 = -t771 * t804 + t774 * t807;
t845 = t747 * t756;
t748 = -t772 * t805 + t775 * t808;
t844 = t748 * t757;
t843 = t755 * t773;
t842 = t756 * t774;
t841 = t757 * t775;
t840 = t803 * t806;
t839 = t804 * t807;
t838 = t805 * t808;
t776 = pkin(2) * t851;
t834 = t813 * t788 + t776;
t777 = pkin(2) * t853;
t833 = t814 * t789 + t777;
t778 = pkin(2) * t855;
t832 = t815 * t790 + t778;
t828 = pkin(2) * t854;
t827 = pkin(2) * t852;
t826 = pkin(2) * t850;
t825 = t787 * t815 - t828;
t824 = t786 * t814 - t827;
t823 = t785 * t813 - t826;
t821 = koppelP(1,1);
t820 = koppelP(2,1);
t819 = koppelP(3,1);
t818 = koppelP(1,2);
t817 = koppelP(2,2);
t816 = koppelP(3,2);
t812 = rSges(4,1);
t811 = rSges(4,2);
t810 = xP(3);
t793 = m(1) + m(2) + m(3);
t792 = cos(t810);
t791 = sin(t810);
t769 = -t791 * t818 + t792 * t821;
t768 = -t791 * t817 + t792 * t820;
t767 = -t791 * t816 + t792 * t819;
t766 = -t791 * t821 - t792 * t818;
t765 = -t791 * t820 - t792 * t817;
t764 = -t791 * t819 - t792 * t816;
t763 = t784 * t787 + 0.2e1 * t828;
t762 = t783 * t786 + 0.2e1 * t827;
t761 = t782 * t785 + 0.2e1 * t826;
t760 = t790 * t784 - 0.2e1 * t778;
t759 = t789 * t783 - 0.2e1 * t777;
t758 = t788 * t782 - 0.2e1 * t776;
t754 = -0.2e1 * t790 * t835 + t805 * (pkin(2) * t790 - t855);
t753 = t787 * t829 - t805 * (pkin(2) * t787 + t854);
t752 = -0.2e1 * t789 * t836 + t804 * (pkin(2) * t789 - t853);
t751 = t786 * t830 - t804 * (pkin(2) * t786 + t852);
t750 = -0.2e1 * t788 * t837 + t803 * (pkin(2) * t788 - t851);
t749 = t785 * t831 - t803 * (pkin(2) * t785 + t850);
t745 = t825 * t808 - t805 * (t790 + t832);
t744 = t824 * t807 - t804 * (t789 + t833);
t743 = t823 * t806 - t803 * (t788 + t834);
t742 = t832 * t808 - t805 * (-t787 - t825);
t741 = t833 * t807 - t804 * (-t786 - t824);
t740 = t834 * t806 - t803 * (-t785 - t823);
t739 = -t760 * t838 + t763 * t796 - t787 * t822 - t787 - t828;
t738 = -t759 * t839 + t762 * t795 - t786 * t822 - t786 - t827;
t737 = -t758 * t840 + t761 * t794 - t785 * t822 - t785 - t826;
t736 = t760 * t796 + t763 * t838 - t822 * t790 + t778 - t790;
t735 = t759 * t795 + t762 * t839 - t822 * t789 + t777 - t789;
t734 = t758 * t794 + t761 * t840 - t822 * t788 + t776 - t788;
t1 = [t749 * t849 + t751 * t848 + t753 * t847 - m(4) * g(1) + (-t734 * t843 - t735 * t842 - t736 * t841) * t793 + (t740 * t846 + t741 * t845 + t742 * t844) * m(3); t750 * t849 + t752 * t848 + t754 * t847 - m(4) * g(2) + (-t737 * t843 - t738 * t842 - t739 * t841) * t793 + (t743 * t846 + t744 * t845 + t745 * t844) * m(3); m(4) * ((g(1) * t812 + g(2) * t811) * t791 + (g(1) * t811 - g(2) * t812) * t792) + (-(t736 * t766 + t739 * t769) * t775 * t793 + (t753 * t766 + t754 * t769) * t733 + (t742 * t766 + t745 * t769) * m(3) * t748) * t757 + (-(t735 * t765 + t738 * t768) * t774 * t793 + (t751 * t765 + t752 * t768) * t732 + (t741 * t765 + t744 * t768) * m(3) * t747) * t756 + (-(t734 * t764 + t737 * t767) * t773 * t793 + (t749 * t764 + t750 * t767) * t731 + (t740 * t764 + t743 * t767) * m(3) * t746) * t755;];
taugX  = t1;
