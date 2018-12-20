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
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3PRP1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:35:05
% EndTime: 2018-12-20 17:35:05
% DurationCPUTime: 0.59s
% Computational Cost: add. (654->143), mult. (1133->231), div. (36->3), fcn. (572->14), ass. (0->126)
t835 = 2 * pkin(2);
t793 = pkin(2) ^ 2;
t834 = 1 + t793;
t833 = mrSges(3,3) - mrSges(2,2);
t774 = legFrame(1,3);
t762 = sin(t774);
t832 = qJ(3,1) * t762;
t765 = cos(t774);
t831 = qJ(3,1) * t765;
t773 = legFrame(2,3);
t761 = sin(t773);
t830 = qJ(3,2) * t761;
t764 = cos(t773);
t829 = qJ(3,2) * t764;
t772 = legFrame(3,3);
t760 = sin(t772);
t828 = qJ(3,3) * t760;
t763 = cos(t772);
t827 = qJ(3,3) * t763;
t747 = -t760 * g(1) + t763 * g(2);
t750 = t763 * g(1) + t760 * g(2);
t753 = m(3) * qJ(3,3) + t833;
t756 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t775 = sin(qJ(2,3));
t778 = cos(qJ(2,3));
t711 = (-t747 * t756 - t753 * t750) * t778 - (t747 * t753 - t756 * t750) * t775;
t784 = qJ(3,3) ^ 2;
t757 = -t784 + t834;
t769 = t778 ^ 2;
t814 = t778 * qJ(3,3);
t735 = 0.1e1 / (t775 * t814 * t835 + t757 * t769 - t784 - t834);
t826 = t711 * t735;
t748 = -t761 * g(1) + t764 * g(2);
t751 = t764 * g(1) + t761 * g(2);
t754 = m(3) * qJ(3,2) + t833;
t776 = sin(qJ(2,2));
t779 = cos(qJ(2,2));
t712 = (-t748 * t756 - t754 * t751) * t779 - (t748 * t754 - t756 * t751) * t776;
t785 = qJ(3,2) ^ 2;
t758 = -t785 + t834;
t770 = t779 ^ 2;
t813 = t779 * qJ(3,2);
t736 = 0.1e1 / (t776 * t813 * t835 + t758 * t770 - t785 - t834);
t825 = t712 * t736;
t749 = -t762 * g(1) + t765 * g(2);
t752 = t765 * g(1) + t762 * g(2);
t755 = m(3) * qJ(3,1) + t833;
t777 = sin(qJ(2,1));
t780 = cos(qJ(2,1));
t713 = (-t749 * t756 - t755 * t752) * t780 - (t749 * t755 - t756 * t752) * t777;
t786 = qJ(3,1) ^ 2;
t759 = -t786 + t834;
t771 = t780 ^ 2;
t812 = t780 * qJ(3,1);
t737 = 0.1e1 / (t777 * t812 * t835 + t759 * t771 - t786 - t834);
t824 = t713 * t737;
t726 = -t778 * t747 + t775 * t750;
t823 = t726 * t735;
t727 = -t779 * t748 + t776 * t751;
t822 = t727 * t736;
t728 = -t780 * t749 + t777 * t752;
t821 = t728 * t737;
t820 = t735 * t747;
t819 = t736 * t748;
t818 = t737 * t749;
t817 = t775 * t778;
t816 = t776 * t779;
t815 = t777 * t780;
t811 = -0.2e1 * t814;
t810 = -0.2e1 * t813;
t809 = -0.2e1 * t812;
t808 = pkin(2) * t832;
t807 = pkin(2) * t830;
t806 = pkin(2) * t828;
t805 = pkin(2) * t827;
t804 = pkin(2) * t829;
t803 = pkin(2) * t831;
t802 = -t762 * t759 + 0.2e1 * t803;
t801 = -t786 * t762 - t803;
t800 = t765 * t786 - t808;
t799 = -t761 * t758 + 0.2e1 * t804;
t798 = -t785 * t761 - t804;
t797 = t764 * t785 - t807;
t796 = -t760 * t757 + 0.2e1 * t805;
t795 = -t784 * t760 - t805;
t794 = t763 * t784 - t806;
t792 = koppelP(1,1);
t791 = koppelP(2,1);
t790 = koppelP(3,1);
t789 = koppelP(1,2);
t788 = koppelP(2,2);
t787 = koppelP(3,2);
t783 = mrSges(4,1);
t782 = mrSges(4,2);
t781 = xP(3);
t768 = m(1) + m(2) + m(3);
t767 = cos(t781);
t766 = sin(t781);
t746 = -t766 * t789 + t767 * t792;
t745 = -t766 * t788 + t767 * t791;
t744 = -t766 * t787 + t767 * t790;
t743 = -t766 * t792 - t767 * t789;
t742 = -t766 * t791 - t767 * t788;
t741 = -t766 * t790 - t767 * t787;
t740 = t759 * t765 + 0.2e1 * t808;
t739 = t758 * t764 + 0.2e1 * t807;
t738 = t757 * t763 + 0.2e1 * t806;
t734 = t765 * t809 + t777 * (pkin(2) * t765 + t832);
t733 = t762 * t809 + t777 * (pkin(2) * t762 - t831);
t732 = t764 * t810 + t776 * (pkin(2) * t764 + t830);
t731 = t761 * t810 + t776 * (pkin(2) * t761 - t829);
t730 = t763 * t811 + t775 * (pkin(2) * t763 + t828);
t729 = t760 * t811 + t775 * (pkin(2) * t760 - t827);
t725 = t800 * t780 - t777 * (t762 - t801);
t724 = t797 * t779 - t776 * (t761 - t798);
t723 = t794 * t778 - t775 * (t760 - t795);
t722 = t801 * t780 + t777 * (-t765 - t800);
t721 = t798 * t779 + t776 * (-t764 - t797);
t720 = t795 * t778 + t775 * (-t763 - t794);
t719 = -t740 * t815 + t793 * t762 + t802 * t771 + t762 - t803;
t718 = t740 * t771 - t765 * t793 + t802 * t815 - t765 - t808;
t717 = -t739 * t816 + t793 * t761 + t799 * t770 + t761 - t804;
t716 = t739 * t770 - t764 * t793 + t799 * t816 - t764 - t807;
t715 = -t738 * t817 + t793 * t760 + t796 * t769 + t760 - t805;
t714 = t738 * t769 - t763 * t793 + t796 * t817 - t763 - t806;
t1 = [t730 * t826 + t732 * t825 + t734 * t824 - g(1) * m(4) + (-t715 * t820 - t717 * t819 - t719 * t818) * t768 + (-t720 * t823 - t721 * t822 - t722 * t821) * m(3); t729 * t826 + t731 * t825 + t733 * t824 - g(2) * m(4) + (-t714 * t820 - t716 * t819 - t718 * t818) * t768 + (-t723 * t823 - t724 * t822 - t725 * t821) * m(3); -(-g(1) * t783 - g(2) * t782) * t766 + t767 * (g(1) * t782 - g(2) * t783) + (-(t718 * t746 + t719 * t743) * t749 * t768 + (t733 * t746 + t734 * t743) * t713 - (t722 * t743 + t725 * t746) * m(3) * t728) * t737 + (-(t716 * t745 + t717 * t742) * t748 * t768 + (t731 * t745 + t732 * t742) * t712 - (t721 * t742 + t724 * t745) * m(3) * t727) * t736 + (-(t714 * t744 + t715 * t741) * t747 * t768 + (t729 * t744 + t730 * t741) * t711 - (t720 * t741 + t723 * t744) * m(3) * t726) * t735;];
taugX  = t1;
