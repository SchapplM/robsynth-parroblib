% Calculate Gravitation load for parallel robot
% P3PRP1A0
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
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3PRP1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:34:54
% EndTime: 2018-12-20 17:34:55
% DurationCPUTime: 0.57s
% Computational Cost: add. (690->153), mult. (1242->252), div. (36->3), fcn. (644->14), ass. (0->125)
t849 = 2 * pkin(2);
t808 = pkin(2) ^ 2;
t848 = 1 + t808;
t788 = legFrame(1,3);
t773 = sin(t788);
t847 = qJ(3,1) * t773;
t776 = cos(t788);
t846 = qJ(3,1) * t776;
t787 = legFrame(2,3);
t772 = sin(t787);
t845 = qJ(3,2) * t772;
t775 = cos(t787);
t844 = qJ(3,2) * t775;
t786 = legFrame(3,3);
t771 = sin(t786);
t843 = qJ(3,3) * t771;
t774 = cos(t786);
t842 = qJ(3,3) * t774;
t762 = -t771 * g(1) + t774 * g(2);
t765 = t774 * g(1) + t771 * g(2);
t783 = rSges(3,3) + qJ(3,3);
t789 = sin(qJ(2,3));
t792 = cos(qJ(2,3));
t795 = pkin(2) + rSges(3,1);
t726 = ((-t762 * t795 - t783 * t765) * m(3) + m(2) * (-rSges(2,1) * t762 + rSges(2,2) * t765)) * t792 + t789 * ((-t762 * t783 + t795 * t765) * m(3) + m(2) * (rSges(2,1) * t765 + rSges(2,2) * t762));
t799 = qJ(3,3) ^ 2;
t768 = -t799 + t848;
t780 = t792 ^ 2;
t829 = t792 * qJ(3,3);
t750 = 0.1e1 / (t789 * t829 * t849 + t768 * t780 - t799 - t848);
t841 = t726 * t750;
t763 = -t772 * g(1) + t775 * g(2);
t766 = t775 * g(1) + t772 * g(2);
t784 = rSges(3,3) + qJ(3,2);
t790 = sin(qJ(2,2));
t793 = cos(qJ(2,2));
t727 = ((-t763 * t795 - t784 * t766) * m(3) + m(2) * (-rSges(2,1) * t763 + rSges(2,2) * t766)) * t793 + t790 * ((-t763 * t784 + t795 * t766) * m(3) + m(2) * (rSges(2,1) * t766 + rSges(2,2) * t763));
t800 = qJ(3,2) ^ 2;
t769 = -t800 + t848;
t781 = t793 ^ 2;
t828 = t793 * qJ(3,2);
t751 = 0.1e1 / (t790 * t828 * t849 + t769 * t781 - t800 - t848);
t840 = t727 * t751;
t764 = -t773 * g(1) + t776 * g(2);
t767 = t776 * g(1) + t773 * g(2);
t785 = rSges(3,3) + qJ(3,1);
t791 = sin(qJ(2,1));
t794 = cos(qJ(2,1));
t728 = ((-t764 * t795 - t785 * t767) * m(3) + m(2) * (-rSges(2,1) * t764 + rSges(2,2) * t767)) * t794 + t791 * ((-t764 * t785 + t795 * t767) * m(3) + m(2) * (rSges(2,1) * t767 + rSges(2,2) * t764));
t801 = qJ(3,1) ^ 2;
t770 = -t801 + t848;
t782 = t794 ^ 2;
t827 = t794 * qJ(3,1);
t752 = 0.1e1 / (t791 * t827 * t849 + t770 * t782 - t801 - t848);
t839 = t728 * t752;
t741 = -t762 * t792 + t765 * t789;
t838 = t741 * t750;
t742 = -t763 * t793 + t766 * t790;
t837 = t742 * t751;
t743 = -t764 * t794 + t767 * t791;
t836 = t743 * t752;
t835 = t750 * t762;
t834 = t751 * t763;
t833 = t752 * t764;
t832 = t789 * t792;
t831 = t790 * t793;
t830 = t791 * t794;
t826 = -0.2e1 * t829;
t825 = -0.2e1 * t828;
t824 = -0.2e1 * t827;
t823 = pkin(2) * t847;
t822 = pkin(2) * t845;
t821 = pkin(2) * t843;
t820 = pkin(2) * t842;
t819 = pkin(2) * t844;
t818 = pkin(2) * t846;
t817 = -t773 * t770 + 0.2e1 * t818;
t816 = -t801 * t773 - t818;
t815 = t776 * t801 - t823;
t814 = -t772 * t769 + 0.2e1 * t819;
t813 = -t800 * t772 - t819;
t812 = t775 * t800 - t822;
t811 = -t771 * t768 + 0.2e1 * t820;
t810 = -t799 * t771 - t820;
t809 = t774 * t799 - t821;
t807 = koppelP(1,1);
t806 = koppelP(2,1);
t805 = koppelP(3,1);
t804 = koppelP(1,2);
t803 = koppelP(2,2);
t802 = koppelP(3,2);
t798 = rSges(4,1);
t797 = rSges(4,2);
t796 = xP(3);
t779 = m(1) + m(2) + m(3);
t778 = cos(t796);
t777 = sin(t796);
t761 = -t777 * t804 + t778 * t807;
t760 = -t777 * t803 + t778 * t806;
t759 = -t777 * t802 + t778 * t805;
t758 = -t777 * t807 - t778 * t804;
t757 = -t777 * t806 - t778 * t803;
t756 = -t777 * t805 - t778 * t802;
t755 = t770 * t776 + 0.2e1 * t823;
t754 = t769 * t775 + 0.2e1 * t822;
t753 = t768 * t774 + 0.2e1 * t821;
t749 = t776 * t824 + t791 * (pkin(2) * t776 + t847);
t748 = t773 * t824 + t791 * (pkin(2) * t773 - t846);
t747 = t775 * t825 + t790 * (pkin(2) * t775 + t845);
t746 = t772 * t825 + t790 * (pkin(2) * t772 - t844);
t745 = t774 * t826 + t789 * (pkin(2) * t774 + t843);
t744 = t771 * t826 + t789 * (pkin(2) * t771 - t842);
t740 = t815 * t794 - t791 * (t773 - t816);
t739 = t812 * t793 - t790 * (t772 - t813);
t738 = t809 * t792 - t789 * (t771 - t810);
t737 = t816 * t794 + t791 * (-t776 - t815);
t736 = t813 * t793 + t790 * (-t775 - t812);
t735 = t810 * t792 + t789 * (-t774 - t809);
t734 = -t755 * t830 + t808 * t773 + t817 * t782 + t773 - t818;
t733 = t755 * t782 - t776 * t808 + t817 * t830 - t776 - t823;
t732 = -t754 * t831 + t808 * t772 + t814 * t781 + t772 - t819;
t731 = t754 * t781 - t775 * t808 + t814 * t831 - t775 - t822;
t730 = -t753 * t832 + t808 * t771 + t811 * t780 + t771 - t820;
t729 = t753 * t780 - t774 * t808 + t811 * t832 - t774 - t821;
t1 = [t745 * t841 + t747 * t840 + t749 * t839 - m(4) * g(1) + (-t730 * t835 - t732 * t834 - t734 * t833) * t779 + (-t735 * t838 - t736 * t837 - t737 * t836) * m(3); t744 * t841 + t746 * t840 + t748 * t839 - m(4) * g(2) + (-t729 * t835 - t731 * t834 - t733 * t833) * t779 + (-t738 * t838 - t739 * t837 - t740 * t836) * m(3); m(4) * ((g(1) * t798 + g(2) * t797) * t777 + (g(1) * t797 - g(2) * t798) * t778) + (-(t733 * t761 + t734 * t758) * t764 * t779 + (t748 * t761 + t749 * t758) * t728 - (t737 * t758 + t740 * t761) * m(3) * t743) * t752 + (-(t731 * t760 + t732 * t757) * t763 * t779 + (t746 * t760 + t747 * t757) * t727 - (t736 * t757 + t739 * t760) * m(3) * t742) * t751 + (-(t729 * t759 + t730 * t756) * t762 * t779 + (t744 * t759 + t745 * t756) * t726 - (t735 * t756 + t738 * t759) * m(3) * t741) * t750;];
taugX  = t1;
