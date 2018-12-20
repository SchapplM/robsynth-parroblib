% Calculate Gravitation load for parallel robot
% P3RPP1A0
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
%   pkin=[a2,a3,d1]';
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
% Datum: 2018-12-20 17:50
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3RPP1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:50:37
% EndTime: 2018-12-20 17:50:38
% DurationCPUTime: 0.62s
% Computational Cost: add. (765->129), mult. (977->210), div. (36->3), fcn. (452->14), ass. (0->131)
t852 = 2 * pkin(1);
t810 = m(2) + m(3);
t851 = pkin(1) ^ 2 + 1;
t798 = legFrame(3,3);
t784 = sin(t798);
t787 = cos(t798);
t760 = -g(1) * t784 + g(2) * t787;
t763 = g(1) * t787 + g(2) * t784;
t835 = mrSges(3,2) + mrSges(2,3) - mrSges(1,2);
t775 = qJ(2,3) * t810 + t835;
t804 = sin(qJ(1,3));
t807 = cos(qJ(1,3));
t801 = -qJ(3,3) - pkin(1);
t831 = m(2) * pkin(1) + mrSges(1,1) - mrSges(2,2) + mrSges(3,3);
t830 = m(3) * t801 - t831;
t715 = (t760 * t830 - t775 * t763) * t807 - t804 * (t760 * t775 + t830 * t763);
t778 = t851 + (t852 + qJ(3,3)) * qJ(3,3);
t816 = (qJ(2,3) ^ 2);
t772 = 1 / (t816 + t778);
t850 = t715 * t772;
t799 = legFrame(2,3);
t785 = sin(t799);
t788 = cos(t799);
t761 = -g(1) * t785 + g(2) * t788;
t764 = g(1) * t788 + g(2) * t785;
t776 = qJ(2,2) * t810 + t835;
t805 = sin(qJ(1,2));
t808 = cos(qJ(1,2));
t802 = -qJ(3,2) - pkin(1);
t829 = m(3) * t802 - t831;
t716 = (t761 * t829 - t776 * t764) * t808 - t805 * (t761 * t776 + t829 * t764);
t779 = t851 + (t852 + qJ(3,2)) * qJ(3,2);
t818 = (qJ(2,2) ^ 2);
t773 = 1 / (t818 + t779);
t849 = t716 * t773;
t800 = legFrame(1,3);
t786 = sin(t800);
t789 = cos(t800);
t762 = -g(1) * t786 + g(2) * t789;
t765 = g(1) * t789 + g(2) * t786;
t777 = qJ(2,1) * t810 + t835;
t806 = sin(qJ(1,1));
t809 = cos(qJ(1,1));
t803 = -qJ(3,1) - pkin(1);
t828 = m(3) * t803 - t831;
t717 = (t762 * t828 - t777 * t765) * t809 - t806 * (t762 * t777 + t828 * t765);
t780 = t851 + (t852 + qJ(3,1)) * qJ(3,1);
t820 = (qJ(2,1) ^ 2);
t774 = 1 / (t820 + t780);
t848 = t717 * t774;
t736 = t760 * t807 - t763 * t804;
t847 = t736 * t772;
t737 = t760 * t804 + t763 * t807;
t846 = t737 * t772;
t738 = t761 * t808 - t764 * t805;
t845 = t738 * t773;
t739 = t761 * t805 + t764 * t808;
t844 = t739 * t773;
t740 = t762 * t809 - t765 * t806;
t843 = t740 * t774;
t741 = t762 * t806 + t765 * t809;
t842 = t741 * t774;
t841 = t801 * t804;
t840 = t801 * t807;
t839 = t802 * t805;
t838 = t802 * t808;
t837 = t803 * t806;
t836 = t803 * t809;
t834 = qJ(2,1) * t836;
t833 = qJ(2,2) * t838;
t832 = qJ(2,3) * t840;
t826 = koppelP(1,1);
t825 = koppelP(2,1);
t824 = koppelP(3,1);
t823 = koppelP(1,2);
t822 = koppelP(2,2);
t821 = koppelP(3,2);
t814 = mrSges(4,1);
t813 = mrSges(4,2);
t812 = xP(3);
t797 = 1 + t820;
t796 = 1 + t818;
t795 = 1 + t816;
t791 = cos(t812);
t790 = sin(t812);
t783 = qJ(2,1) * t837;
t782 = qJ(2,2) * t839;
t781 = qJ(2,3) * t841;
t771 = qJ(2,1) * t806 - t836;
t770 = qJ(2,2) * t805 - t838;
t769 = qJ(2,3) * t804 - t840;
t768 = qJ(2,1) * t809 + t837;
t767 = qJ(2,2) * t808 + t839;
t766 = qJ(2,3) * t807 + t841;
t759 = -t790 * t823 + t791 * t826;
t758 = -t790 * t822 + t791 * t825;
t757 = -t790 * t821 + t791 * t824;
t756 = -t790 * t826 - t791 * t823;
t755 = -t790 * t825 - t791 * t822;
t754 = -t790 * t824 - t791 * t821;
t753 = -t797 * t809 - t783;
t752 = -t796 * t808 - t782;
t751 = -t795 * t807 - t781;
t750 = t797 * t806 - t834;
t749 = t796 * t805 - t833;
t748 = t795 * t804 - t832;
t747 = t780 * t809 - t783;
t746 = t779 * t808 - t782;
t745 = t778 * t807 - t781;
t744 = t780 * t806 + t834;
t743 = t779 * t805 + t833;
t742 = t778 * t804 + t832;
t735 = t768 * t786 + t771 * t789;
t734 = t767 * t785 + t770 * t788;
t733 = t766 * t784 + t769 * t787;
t732 = t768 * t789 - t771 * t786;
t731 = t767 * t788 - t770 * t785;
t730 = t766 * t787 - t769 * t784;
t729 = t750 * t786 + t753 * t789;
t728 = t749 * t785 + t752 * t788;
t727 = t748 * t784 + t751 * t787;
t726 = t750 * t789 - t753 * t786;
t725 = t749 * t788 - t752 * t785;
t724 = t748 * t787 - t751 * t784;
t723 = -t744 * t786 + t747 * t789;
t722 = -t743 * t785 + t746 * t788;
t721 = -t742 * t784 + t745 * t787;
t720 = t744 * t789 + t747 * t786;
t719 = t743 * t788 + t746 * t785;
t718 = t742 * t787 + t745 * t784;
t1 = [t730 * t850 + t731 * t849 + t732 * t848 - g(1) * m(4) + (t724 * t847 + t725 * t845 + t726 * t843) * t810 + (-t721 * t846 - t722 * t844 - t723 * t842) * m(3); t733 * t850 + t734 * t849 + t735 * t848 - g(2) * m(4) + (t727 * t847 + t728 * t845 + t729 * t843) * t810 + (-t718 * t846 - t719 * t844 - t720 * t842) * m(3); -(-g(1) * t814 - g(2) * t813) * t790 + t791 * (g(1) * t813 - g(2) * t814) + ((t732 * t756 + t735 * t759) * t717 + (t726 * t756 + t729 * t759) * t740 * t810 - (t720 * t759 + t723 * t756) * m(3) * t741) * t774 + ((t731 * t755 + t734 * t758) * t716 + (t725 * t755 + t728 * t758) * t738 * t810 - (t719 * t758 + t722 * t755) * m(3) * t739) * t773 + ((t730 * t754 + t733 * t757) * t715 + (t724 * t754 + t727 * t757) * t736 * t810 - (t718 * t757 + t721 * t754) * m(3) * t737) * t772;];
taugX  = t1;
