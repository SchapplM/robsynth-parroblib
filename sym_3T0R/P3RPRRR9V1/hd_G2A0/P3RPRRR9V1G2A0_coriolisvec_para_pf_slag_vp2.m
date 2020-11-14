% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR9V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:39
% EndTime: 2020-08-06 18:51:45
% DurationCPUTime: 5.25s
% Computational Cost: add. (20541->458), mult. (26184->743), div. (4311->10), fcn. (22560->36), ass. (0->303)
t1006 = -2 * pkin(5);
t822 = sin(pkin(7));
t1000 = -0.2e1 * t822;
t823 = cos(pkin(7));
t781 = t823 * pkin(2);
t1005 = -0.2e1 * t781;
t1003 = 2 * pkin(1);
t1004 = t1003 / 0.2e1;
t768 = t781 + pkin(1);
t813 = pkin(7) + qJ(3,1);
t780 = cos(t813);
t992 = pkin(3) * t780;
t746 = t768 + t992;
t812 = pkin(7) + qJ(3,2);
t779 = cos(t812);
t993 = pkin(3) * t779;
t745 = t768 + t993;
t811 = pkin(7) + qJ(3,3);
t778 = cos(t811);
t994 = pkin(3) * t778;
t744 = t768 + t994;
t824 = pkin(5) + qJ(2,3);
t762 = mrSges(3,1) * t824 - Ifges(3,5);
t825 = pkin(5) + qJ(2,2);
t763 = mrSges(3,1) * t825 - Ifges(3,5);
t826 = pkin(5) + qJ(2,1);
t764 = mrSges(3,1) * t826 - Ifges(3,5);
t759 = mrSges(3,2) * t824 - Ifges(3,6);
t760 = mrSges(3,2) * t825 - Ifges(3,6);
t761 = mrSges(3,2) * t826 - Ifges(3,6);
t1002 = 0.2e1 * pkin(3);
t810 = t823 ^ 2;
t1001 = 0.2e1 * t810;
t999 = 0.4e1 * t823;
t998 = -4 * pkin(5) - 4 * pkin(6);
t997 = pkin(1) * mrSges(3,2);
t996 = pkin(2) * mrSges(3,1);
t995 = pkin(2) * mrSges(3,2);
t775 = sin(t811);
t829 = legFrame(3,2);
t795 = cos(t829);
t792 = sin(t829);
t833 = sin(qJ(1,3));
t955 = t792 * t833;
t732 = t775 * t795 - t778 * t955;
t952 = t795 * t833;
t733 = t775 * t792 + t778 * t952;
t772 = 0.1e1 / t778;
t807 = pkin(6) + t824;
t789 = 1 / t807;
t839 = cos(qJ(1,3));
t845 = xDP(3);
t846 = xDP(2);
t847 = xDP(1);
t705 = (t839 * t845 + (t732 * t846 + t733 * t847) * t772) * t789;
t703 = pkin(1) * t705;
t776 = sin(t812);
t830 = legFrame(2,2);
t796 = cos(t830);
t793 = sin(t830);
t835 = sin(qJ(1,2));
t954 = t793 * t835;
t734 = t776 * t796 - t779 * t954;
t951 = t796 * t835;
t735 = t776 * t793 + t779 * t951;
t773 = 0.1e1 / t779;
t808 = pkin(6) + t825;
t790 = 1 / t808;
t841 = cos(qJ(1,2));
t706 = (t841 * t845 + (t734 * t846 + t735 * t847) * t773) * t790;
t704 = pkin(1) * t706;
t782 = pkin(1) * t822;
t777 = sin(t813);
t831 = legFrame(1,2);
t797 = cos(t831);
t794 = sin(t831);
t837 = sin(qJ(1,1));
t953 = t794 * t837;
t736 = t777 * t797 - t780 * t953;
t950 = t797 * t837;
t737 = t777 * t794 + t780 * t950;
t774 = 0.1e1 / t780;
t809 = pkin(6) + t826;
t791 = 1 / t809;
t843 = cos(qJ(1,1));
t707 = (t843 * t845 + (t736 * t846 + t737 * t847) * t774) * t791;
t702 = t707 * pkin(1);
t838 = cos(qJ(3,3));
t818 = t838 ^ 2;
t991 = t818 * pkin(3);
t840 = cos(qJ(3,2));
t819 = t840 ^ 2;
t990 = t819 * pkin(3);
t842 = cos(qJ(3,1));
t820 = t842 ^ 2;
t989 = t820 * pkin(3);
t988 = t838 * pkin(2);
t987 = t840 * pkin(2);
t986 = t842 * pkin(2);
t828 = mrSges(2,3) + mrSges(3,3);
t827 = Ifges(3,1) - Ifges(3,2);
t985 = Ifges(2,4) - Ifges(3,4);
t984 = mrSges(3,1) * t822;
t890 = pkin(1) * t837 - t809 * t843;
t836 = sin(qJ(3,1));
t939 = t836 * t822;
t719 = t890 * t939 + (t820 - 0.1e1) * t837 * pkin(3);
t731 = pkin(1) * t836 + (-pkin(3) + t986 + 0.2e1 * t989) * t822;
t862 = -pkin(3) / 0.2e1;
t749 = t989 + t986 / 0.2e1 + t862;
t899 = t837 * t939;
t872 = pkin(2) * t899 + (t899 * t1002 - t890) * t842;
t935 = t842 * (-pkin(3) * t836 + t782);
t863 = pkin(2) / 0.2e1;
t956 = (pkin(3) * t842 + t863) * t836;
t697 = (-t749 * t953 + t797 * t956) * t1001 + (t797 * t731 + t794 * t872) * t823 + t719 * t794 + t797 * t935;
t743 = 0.1e1 / (t842 * t823 - t939);
t959 = t743 * t791;
t688 = t697 * t846 * t959;
t698 = (t749 * t950 + t794 * t956) * t1001 + (t794 * t731 - t797 * t872) * t823 - t719 * t797 + t794 * t935;
t689 = t698 * t847 * t959;
t728 = t746 * t843 + t809 * t837;
t716 = t728 * t791 * t845;
t678 = t702 - 0.2e1 * t689 - 0.2e1 * t688 - 0.2e1 * t716;
t983 = mrSges(3,2) * t678;
t892 = pkin(1) * t833 - t807 * t839;
t832 = sin(qJ(3,3));
t943 = t832 * t822;
t717 = t892 * t943 + (t818 - 0.1e1) * t833 * pkin(3);
t729 = pkin(1) * t832 + (-pkin(3) + t988 + 0.2e1 * t991) * t822;
t747 = t991 + t988 / 0.2e1 + t862;
t901 = t833 * t943;
t874 = pkin(2) * t901 + (t901 * t1002 - t892) * t838;
t937 = t838 * (-pkin(3) * t832 + t782);
t958 = (pkin(3) * t838 + t863) * t832;
t693 = (-t747 * t955 + t795 * t958) * t1001 + (t795 * t729 + t792 * t874) * t823 + t717 * t792 + t795 * t937;
t741 = 0.1e1 / (t838 * t823 - t943);
t961 = t741 * t789;
t684 = t693 * t846 * t961;
t694 = (t747 * t952 + t792 * t958) * t1001 + (t792 * t729 - t795 * t874) * t823 - t717 * t795 + t792 * t937;
t685 = t694 * t847 * t961;
t726 = t744 * t839 + t807 * t833;
t714 = t726 * t789 * t845;
t679 = t703 - 0.2e1 * t685 - 0.2e1 * t684 - 0.2e1 * t714;
t982 = mrSges(3,2) * t679;
t891 = pkin(1) * t835 - t808 * t841;
t834 = sin(qJ(3,2));
t941 = t834 * t822;
t718 = t891 * t941 + (t819 - 0.1e1) * t835 * pkin(3);
t730 = pkin(1) * t834 + (-pkin(3) + t987 + 0.2e1 * t990) * t822;
t748 = t990 + t987 / 0.2e1 + t862;
t900 = t835 * t941;
t873 = pkin(2) * t900 + (t900 * t1002 - t891) * t840;
t936 = t840 * (-pkin(3) * t834 + t782);
t957 = (pkin(3) * t840 + t863) * t834;
t695 = (-t748 * t954 + t796 * t957) * t1001 + (t796 * t730 + t793 * t873) * t823 + t718 * t793 + t796 * t936;
t742 = 0.1e1 / (t840 * t823 - t941);
t960 = t742 * t790;
t686 = t695 * t846 * t960;
t696 = (t748 * t951 + t793 * t957) * t1001 + (t793 * t730 - t796 * t873) * t823 - t718 * t796 + t793 * t936;
t687 = t696 * t847 * t960;
t727 = t745 * t841 + t808 * t835;
t715 = t727 * t790 * t845;
t680 = t704 - 0.2e1 * t687 - 0.2e1 * t686 - 0.2e1 * t715;
t981 = mrSges(3,2) * t680;
t980 = mrSges(3,2) * t832;
t979 = mrSges(3,2) * t834;
t978 = mrSges(3,2) * t836;
t786 = mrSges(3,2) * t838;
t787 = mrSges(3,2) * t840;
t788 = mrSges(3,2) * t842;
t977 = Ifges(3,4) * t832;
t976 = Ifges(3,4) * t834;
t975 = Ifges(3,4) * t836;
t974 = t818 * Ifges(3,4);
t973 = t819 * Ifges(3,4);
t972 = t820 * Ifges(3,4);
t783 = t832 * mrSges(3,1);
t784 = t834 * mrSges(3,1);
t785 = t836 * mrSges(3,1);
t971 = m(3) * pkin(2) + mrSges(2,1);
t970 = t705 * t810;
t969 = t705 * t822;
t968 = t706 * t810;
t967 = t706 * t822;
t966 = t707 * t810;
t965 = t707 * t822;
t868 = 0.1e1 / pkin(3);
t723 = (t792 * t847 + t795 * t846) * t868 * t772;
t964 = t723 ^ 2 * t772;
t724 = (t793 * t847 + t796 * t846) * t868 * t773;
t963 = t724 ^ 2 * t773;
t725 = (t794 * t847 + t797 * t846) * t868 * t774;
t962 = t725 ^ 2 * t774;
t949 = t827 * t818;
t948 = t827 * t819;
t947 = t827 * t820;
t946 = t827 * t832;
t945 = t827 * t834;
t944 = t827 * t836;
t942 = t832 * t838;
t940 = t834 * t840;
t938 = t836 * t842;
t676 = t703 - t685 / 0.2e1 - t684 / 0.2e1 - t714 / 0.2e1;
t681 = t685 + t684 + t714;
t769 = 0.2e1 * t811;
t859 = 0.2e1 * pkin(7);
t864 = qJ(2,3) ^ 2;
t867 = pkin(3) ^ 2;
t869 = pkin(2) ^ 2;
t870 = pkin(1) ^ 2;
t871 = -t869 * cos(t859) - (2 * pkin(6) ^ 2) - t867 - t869 - (2 * t870) + ((-4 * pkin(6) + t1006) * pkin(5));
t663 = ((t676 * t1005 + ((qJ(2,3) * t998) - t867 * cos(t769) - (2 * t864) + t871) * t705 / 0.2e1 + (t744 + t1004) * t681) * t705 + ((t723 * t807 * t775 - 0.2e1 * t676 * t778 + (-cos(qJ(3,3) + t859) * pkin(2) - t988) * t705) * t705 - (-t705 * t807 * sin(t769) / 0.2e1 + t723 * t744) * t772 * t723) * pkin(3)) * t789;
t672 = (-pkin(3) * t964 + (-t703 + 0.2e1 * t681 + (-t781 - t994) * t705) * t705) * t789;
t699 = t705 * t946;
t708 = -(-t759 * t838 - t762 * t832) * t823 - t822 * (t759 * t832 - t762 * t838);
t849 = m(2) + m(3);
t805 = pkin(1) * t849;
t884 = mrSges(3,1) * t838 - t980;
t931 = t786 + t783;
t711 = -(-t884 - t971) * t823 - (mrSges(2,2) + t931) * t822 + t805;
t911 = m(3) * pkin(5) + t828;
t750 = qJ(2,3) * t849 + t911;
t803 = -t995 / 0.2e1;
t804 = t996 / 0.4e1;
t815 = 0.2e1 * t996;
t850 = Ifges(3,6) / 0.2e1;
t851 = Ifges(3,5) / 0.4e1;
t860 = -pkin(5) / 0.4e1;
t861 = -pkin(5) / 0.2e1;
t875 = (mrSges(3,3) * t1006) - t849 * t870 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3);
t885 = t869 * m(3) - Ifges(2,1) + Ifges(2,2) + t827;
t888 = t949 * t969;
t889 = (mrSges(3,1) * t1003) + t995 * t1000;
t904 = t775 * t964;
t917 = pkin(2) * t980;
t907 = (t917 - t827) * t969;
t910 = t705 * t974;
t913 = -t996 / 0.2e1;
t914 = -t997 / 0.2e1;
t920 = mrSges(3,1) * t703;
t924 = 0.4e1 * t977;
t925 = 0.4e1 * Ifges(3,4) * t822;
t927 = 0.2e1 * t782;
t928 = -2 * t997;
t934 = (-((t815 + t924) * t838 - 0.2e1 * t917 + t885) * t810 - (t818 * t925 + t889 * t838 + (t971 - t980) * t1003) * t823 - t949 + (t783 + mrSges(2,2)) * t927 - (t824 ^ 2) * m(3) - m(2) * t864 + t875) * t672 + t711 * t663 - t708 * t904 + 0.2e1 * (Ifges(3,4) * t942 + t810 * t949 - t828 * qJ(2,3) + (-(-pkin(2) * t783 + t827 * t942 + t985) * t823 + pkin(1) * t786) * t822) * t672 + 0.4e1 * (0.2e1 * t974 + (t803 + t946) * t838 + t832 * t913 - Ifges(3,4)) * t723 * t970 + (t888 + (((-qJ(2,3) / 0.4e1 + t860) * mrSges(3,1) + t851) * t723 + ((t804 + t977) * t1000 + t914) * t705) * t838 + t907 / 0.2e1 - t832 * (-t759 * t723 + 0.2e1 * t920) / 0.4e1) * t723 * t999 - 0.4e1 * t723 * t910 - 0.2e1 * ((((-qJ(2,3) / 0.2e1 + t861) * mrSges(3,2) + t850) * t723 + t920) * t822 + t699) * t723 * t838 - t723 * (t705 * t928 - t762 * t723) * t943 + 0.2e1 * t705 * (Ifges(3,4) * t723 + t681 * t750);
t677 = t704 - t687 / 0.2e1 - t686 / 0.2e1 - t715 / 0.2e1;
t682 = t687 + t686 + t715;
t770 = 0.2e1 * t812;
t865 = qJ(2,2) ^ 2;
t664 = ((t677 * t1005 + ((qJ(2,2) * t998) - t867 * cos(t770) - (2 * t865) + t871) * t706 / 0.2e1 + (t745 + t1004) * t682) * t706 + ((t724 * t808 * t776 - 0.2e1 * t677 * t779 + (-cos(t859 + qJ(3,2)) * pkin(2) - t987) * t706) * t706 - (-t706 * t808 * sin(t770) / 0.2e1 + t724 * t745) * t773 * t724) * pkin(3)) * t790;
t673 = (-pkin(3) * t963 + (-t704 + 0.2e1 * t682 + (-t781 - t993) * t706) * t706) * t790;
t700 = t706 * t945;
t709 = -(-t760 * t840 - t763 * t834) * t823 - t822 * (t760 * t834 - t763 * t840);
t883 = mrSges(3,1) * t840 - t979;
t930 = t787 + t784;
t712 = -(-t883 - t971) * t823 - (mrSges(2,2) + t930) * t822 + t805;
t751 = qJ(2,2) * t849 + t911;
t887 = t948 * t967;
t903 = t776 * t963;
t916 = pkin(2) * t979;
t906 = (t916 - t827) * t967;
t909 = t706 * t973;
t919 = mrSges(3,1) * t704;
t923 = 0.4e1 * t976;
t933 = (-((t815 + t923) * t840 - 0.2e1 * t916 + t885) * t810 - (t819 * t925 + t889 * t840 + (t971 - t979) * t1003) * t823 - t948 + (t784 + mrSges(2,2)) * t927 - (t825 ^ 2) * m(3) - m(2) * t865 + t875) * t673 + t712 * t664 - t709 * t903 + 0.2e1 * (Ifges(3,4) * t940 + t810 * t948 - t828 * qJ(2,2) + (-(-pkin(2) * t784 + t827 * t940 + t985) * t823 + pkin(1) * t787) * t822) * t673 + 0.4e1 * (0.2e1 * t973 + (t803 + t945) * t840 + t834 * t913 - Ifges(3,4)) * t724 * t968 + (t887 + (((-qJ(2,2) / 0.4e1 + t860) * mrSges(3,1) + t851) * t724 + ((t804 + t976) * t1000 + t914) * t706) * t840 + t906 / 0.2e1 - t834 * (-t760 * t724 + 0.2e1 * t919) / 0.4e1) * t724 * t999 - 0.4e1 * t724 * t909 - 0.2e1 * ((((-qJ(2,2) / 0.2e1 + t861) * mrSges(3,2) + t850) * t724 + t919) * t822 + t700) * t724 * t840 - t724 * (t706 * t928 - t763 * t724) * t941 + 0.2e1 * t706 * (Ifges(3,4) * t724 + t682 * t751);
t675 = t702 - t689 / 0.2e1 - t688 / 0.2e1 - t716 / 0.2e1;
t683 = t689 + t688 + t716;
t771 = 0.2e1 * t813;
t866 = qJ(2,1) ^ 2;
t665 = ((t675 * t1005 + ((qJ(2,1) * t998) - t867 * cos(t771) - (2 * t866) + t871) * t707 / 0.2e1 + (t746 + t1004) * t683) * t707 + ((t725 * t809 * t777 - 0.2e1 * t675 * t780 + (-cos(qJ(3,1) + t859) * pkin(2) - t986) * t707) * t707 - (-t707 * t809 * sin(t771) / 0.2e1 + t725 * t746) * t774 * t725) * pkin(3)) * t791;
t674 = (-pkin(3) * t962 + (-t702 + 0.2e1 * t683 + (-t781 - t992) * t707) * t707) * t791;
t701 = t707 * t944;
t710 = -(-t761 * t842 - t764 * t836) * t823 - t822 * (t761 * t836 - t764 * t842);
t882 = mrSges(3,1) * t842 - t978;
t929 = t788 + t785;
t713 = -(-t882 - t971) * t823 - (mrSges(2,2) + t929) * t822 + t805;
t752 = qJ(2,1) * t849 + t911;
t886 = t947 * t965;
t902 = t777 * t962;
t915 = pkin(2) * t978;
t905 = (t915 - t827) * t965;
t908 = t707 * t972;
t918 = mrSges(3,1) * t702;
t922 = 0.4e1 * t975;
t932 = (-((t815 + t922) * t842 - 0.2e1 * t915 + t885) * t810 - (t820 * t925 + t889 * t842 + (t971 - t978) * t1003) * t823 - t947 + (t785 + mrSges(2,2)) * t927 - (t826 ^ 2) * m(3) - m(2) * t866 + t875) * t674 + t713 * t665 - t710 * t902 + 0.2e1 * (Ifges(3,4) * t938 + t810 * t947 - t828 * qJ(2,1) + (-(-pkin(2) * t785 + t827 * t938 + t985) * t823 + pkin(1) * t788) * t822) * t674 + 0.4e1 * (0.2e1 * t972 + (t803 + t944) * t842 + t836 * t913 - Ifges(3,4)) * t725 * t966 + (t886 + (((-qJ(2,1) / 0.4e1 + t860) * mrSges(3,1) + t851) * t725 + ((t804 + t975) * t1000 + t914) * t707) * t842 + t905 / 0.2e1 - t836 * (-t761 * t725 + 0.2e1 * t918) / 0.4e1) * t725 * t999 - 0.4e1 * t725 * t908 - 0.2e1 * ((((-qJ(2,1) / 0.2e1 + t861) * mrSges(3,2) + t850) * t725 + t918) * t822 + t701) * t725 * t842 - t725 * (t707 * t928 - t764 * t725) * t939 + 0.2e1 * t707 * (Ifges(3,4) * t725 + t683 * t752);
t912 = -t996 / 0.4e1;
t898 = t934 * t772;
t897 = t933 * t773;
t896 = t932 * t774;
t895 = -t663 * t849 + t672 * t711 - (t705 * t750 + 0.2e1 * (-t884 * t822 - t931 * t823) * t723) * t705;
t894 = -t664 * t849 + t673 * t712 - (t706 * t751 + 0.2e1 * (-t883 * t822 - t930 * t823) * t724) * t706;
t893 = -t665 * t849 + t674 * t713 - (t707 * t752 + 0.2e1 * (-t882 * t822 - t929 * t823) * t725) * t707;
t798 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t802 = -t995 / 0.4e1;
t852 = -Ifges(3,4) / 0.2e1;
t881 = t772 * ((-0.4e1 * (t974 + (t798 * t832 + t802) * t838 + t832 * t912 + t852) * t970 + (-0.2e1 * t888 + (t982 + (t924 + t996) * t969) * t838 - t907 + t679 * t783) * t823 + 0.2e1 * t910 + (t679 * t984 + t699) * t838 - t943 * t982 - Ifges(3,4) * t705) * t705 + Ifges(3,3) * t904 + t672 * t708);
t880 = t773 * ((-0.4e1 * (t973 + (t798 * t834 + t802) * t840 + t834 * t912 + t852) * t968 + (-0.2e1 * t887 + (t981 + (t923 + t996) * t967) * t840 - t906 + t680 * t784) * t823 + 0.2e1 * t909 + (t680 * t984 + t700) * t840 - t941 * t981 - Ifges(3,4) * t706) * t706 + Ifges(3,3) * t903 + t673 * t709);
t879 = t774 * ((-0.4e1 * (t972 + (t798 * t836 + t802) * t842 + t836 * t912 + t852) * t966 + (-0.2e1 * t886 + (t983 + (t922 + t996) * t965) * t842 - t905 + t678 * t785) * t823 + 0.2e1 * t908 + (t678 * t984 + t701) * t842 - t939 * t983 - Ifges(3,4) * t707) * t707 + Ifges(3,3) * t902 + t674 * t710);
t878 = t895 * t741;
t877 = t894 * t742;
t876 = t893 * t743;
t1 = [(t698 * t876 + t737 * t896) * t791 + (t696 * t877 + t735 * t897) * t790 + (t694 * t878 + t733 * t898) * t789 + (t792 * t881 + t793 * t880 + t794 * t879) * t868; (t697 * t876 + t736 * t896) * t791 + (t695 * t877 + t734 * t897) * t790 + (t693 * t878 + t732 * t898) * t789 + (t795 * t881 + t796 * t880 + t797 * t879) * t868; (t893 * t728 + t932 * t843) * t791 + (t894 * t727 + t933 * t841) * t790 + (t895 * t726 + t934 * t839) * t789;];
taucX  = t1;
