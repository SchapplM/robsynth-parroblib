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
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-05-03 14:52:15
% EndTime: 2019-05-03 14:52:16
% DurationCPUTime: 0.64s
% Computational Cost: add. (765->129), mult. (977->210), div. (36->3), fcn. (452->14), ass. (0->131)
t854 = 2 * pkin(1);
t812 = m(2) + m(3);
t853 = pkin(1) ^ 2 + 1;
t800 = legFrame(3,3);
t786 = sin(t800);
t789 = cos(t800);
t762 = -g(1) * t786 + g(2) * t789;
t765 = g(1) * t789 + g(2) * t786;
t837 = mrSges(3,2) + mrSges(2,3) - mrSges(1,2);
t777 = qJ(2,3) * t812 + t837;
t806 = sin(qJ(1,3));
t809 = cos(qJ(1,3));
t803 = -qJ(3,3) - pkin(1);
t833 = m(2) * pkin(1) + mrSges(1,1) - mrSges(2,2) + mrSges(3,3);
t832 = m(3) * t803 - t833;
t717 = (t832 * t762 - t777 * t765) * t809 - t806 * (t762 * t777 + t832 * t765);
t780 = t853 + (t854 + qJ(3,3)) * qJ(3,3);
t818 = (qJ(2,3) ^ 2);
t774 = 1 / (t818 + t780);
t852 = t717 * t774;
t801 = legFrame(2,3);
t787 = sin(t801);
t790 = cos(t801);
t763 = -g(1) * t787 + g(2) * t790;
t766 = g(1) * t790 + g(2) * t787;
t778 = qJ(2,2) * t812 + t837;
t807 = sin(qJ(1,2));
t810 = cos(qJ(1,2));
t804 = -qJ(3,2) - pkin(1);
t831 = m(3) * t804 - t833;
t718 = (t831 * t763 - t778 * t766) * t810 - t807 * (t763 * t778 + t831 * t766);
t781 = t853 + (t854 + qJ(3,2)) * qJ(3,2);
t820 = (qJ(2,2) ^ 2);
t775 = 1 / (t820 + t781);
t851 = t718 * t775;
t802 = legFrame(1,3);
t788 = sin(t802);
t791 = cos(t802);
t764 = -g(1) * t788 + g(2) * t791;
t767 = g(1) * t791 + g(2) * t788;
t779 = qJ(2,1) * t812 + t837;
t808 = sin(qJ(1,1));
t811 = cos(qJ(1,1));
t805 = -qJ(3,1) - pkin(1);
t830 = m(3) * t805 - t833;
t719 = (t830 * t764 - t779 * t767) * t811 - t808 * (t764 * t779 + t830 * t767);
t782 = t853 + (t854 + qJ(3,1)) * qJ(3,1);
t822 = (qJ(2,1) ^ 2);
t776 = 1 / (t822 + t782);
t850 = t719 * t776;
t738 = t762 * t809 - t765 * t806;
t849 = t738 * t774;
t739 = t762 * t806 + t765 * t809;
t848 = t739 * t774;
t740 = t763 * t810 - t766 * t807;
t847 = t740 * t775;
t741 = t763 * t807 + t766 * t810;
t846 = t741 * t775;
t742 = t764 * t811 - t767 * t808;
t845 = t742 * t776;
t743 = t764 * t808 + t767 * t811;
t844 = t743 * t776;
t843 = t803 * t806;
t842 = t803 * t809;
t841 = t804 * t807;
t840 = t804 * t810;
t839 = t805 * t808;
t838 = t805 * t811;
t836 = qJ(2,1) * t838;
t835 = qJ(2,2) * t840;
t834 = qJ(2,3) * t842;
t828 = koppelP(1,1);
t827 = koppelP(2,1);
t826 = koppelP(3,1);
t825 = koppelP(1,2);
t824 = koppelP(2,2);
t823 = koppelP(3,2);
t816 = mrSges(4,1);
t815 = mrSges(4,2);
t814 = xP(3);
t799 = 1 + t822;
t798 = 1 + t820;
t797 = 1 + t818;
t793 = cos(t814);
t792 = sin(t814);
t785 = qJ(2,1) * t839;
t784 = qJ(2,2) * t841;
t783 = qJ(2,3) * t843;
t773 = qJ(2,1) * t808 - t838;
t772 = qJ(2,2) * t807 - t840;
t771 = qJ(2,3) * t806 - t842;
t770 = qJ(2,1) * t811 + t839;
t769 = qJ(2,2) * t810 + t841;
t768 = qJ(2,3) * t809 + t843;
t761 = -t792 * t825 + t793 * t828;
t760 = -t792 * t824 + t793 * t827;
t759 = -t792 * t823 + t793 * t826;
t758 = -t792 * t828 - t793 * t825;
t757 = -t792 * t827 - t793 * t824;
t756 = -t792 * t826 - t793 * t823;
t755 = -t799 * t811 - t785;
t754 = -t798 * t810 - t784;
t753 = -t797 * t809 - t783;
t752 = t799 * t808 - t836;
t751 = t798 * t807 - t835;
t750 = t797 * t806 - t834;
t749 = t782 * t811 - t785;
t748 = t781 * t810 - t784;
t747 = t780 * t809 - t783;
t746 = t782 * t808 + t836;
t745 = t781 * t807 + t835;
t744 = t780 * t806 + t834;
t737 = t770 * t788 + t773 * t791;
t736 = t769 * t787 + t772 * t790;
t735 = t768 * t786 + t771 * t789;
t734 = t770 * t791 - t773 * t788;
t733 = t769 * t790 - t772 * t787;
t732 = t768 * t789 - t771 * t786;
t731 = t752 * t788 + t755 * t791;
t730 = t751 * t787 + t754 * t790;
t729 = t750 * t786 + t753 * t789;
t728 = t752 * t791 - t755 * t788;
t727 = t751 * t790 - t754 * t787;
t726 = t750 * t789 - t753 * t786;
t725 = -t746 * t788 + t749 * t791;
t724 = -t745 * t787 + t748 * t790;
t723 = -t744 * t786 + t747 * t789;
t722 = t746 * t791 + t749 * t788;
t721 = t745 * t790 + t748 * t787;
t720 = t744 * t789 + t747 * t786;
t1 = [t732 * t852 + t733 * t851 + t734 * t850 - g(1) * m(4) + (t726 * t849 + t727 * t847 + t728 * t845) * t812 + (-t723 * t848 - t724 * t846 - t725 * t844) * m(3); t735 * t852 + t736 * t851 + t737 * t850 - g(2) * m(4) + (t729 * t849 + t730 * t847 + t731 * t845) * t812 + (-t720 * t848 - t721 * t846 - t722 * t844) * m(3); -(-g(1) * t816 - g(2) * t815) * t792 + t793 * (g(1) * t815 - g(2) * t816) + ((t734 * t758 + t737 * t761) * t719 + (t728 * t758 + t731 * t761) * t742 * t812 - (t722 * t761 + t725 * t758) * m(3) * t743) * t776 + ((t733 * t757 + t736 * t760) * t718 + (t727 * t757 + t730 * t760) * t740 * t812 - (t721 * t760 + t724 * t757) * m(3) * t741) * t775 + ((t732 * t756 + t735 * t759) * t717 + (t726 * t756 + t729 * t759) * t738 * t812 - (t720 * t759 + t723 * t756) * m(3) * t739) * t774;];
taugX  = t1;
