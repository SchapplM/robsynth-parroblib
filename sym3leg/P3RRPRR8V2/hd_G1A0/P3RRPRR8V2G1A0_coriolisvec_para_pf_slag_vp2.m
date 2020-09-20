% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:52
% EndTime: 2020-08-06 21:02:59
% DurationCPUTime: 6.89s
% Computational Cost: add. (27639->455), mult. (34032->702), div. (5265->15), fcn. (25611->62), ass. (0->324)
t853 = -qJ(3,1) - pkin(5);
t1049 = mrSges(3,2) * t853 + Ifges(3,6);
t769 = mrSges(3,1) * t853 + Ifges(3,5);
t845 = sin(pkin(7));
t846 = cos(pkin(7));
t1069 = -t1049 * t845 + t769 * t846;
t852 = -qJ(3,2) - pkin(5);
t1050 = mrSges(3,2) * t852 + Ifges(3,6);
t768 = mrSges(3,1) * t852 + Ifges(3,5);
t1068 = -t1050 * t845 + t768 * t846;
t851 = -qJ(3,3) - pkin(5);
t1051 = mrSges(3,2) * t851 + Ifges(3,6);
t767 = mrSges(3,1) * t851 + Ifges(3,5);
t1067 = -t1051 * t845 + t767 * t846;
t827 = t846 ^ 2;
t1019 = Ifges(3,4) * t827;
t857 = Ifges(3,2) - Ifges(3,1);
t755 = -pkin(2) * mrSges(3,2) - t845 * t857;
t1021 = mrSges(3,1) * t845;
t757 = -pkin(2) * t1021 + Ifges(2,4) - Ifges(3,4);
t715 = t755 * t846 + 0.2e1 * t1019 + t757;
t873 = xDP(3);
t904 = t715 * t873;
t1024 = pkin(3) * t846;
t774 = pkin(2) + t1024;
t866 = cos(qJ(2,3));
t1066 = t774 * t866;
t868 = cos(qJ(2,2));
t1065 = t774 * t868;
t870 = cos(qJ(2,1));
t1064 = t774 * t870;
t833 = qJ(2,1) + pkin(7);
t801 = cos(t833);
t1026 = pkin(3) * t801;
t819 = t870 * pkin(2);
t780 = t819 + pkin(1);
t823 = -pkin(6) + t853;
t865 = sin(qJ(1,1));
t871 = cos(qJ(1,1));
t724 = t780 * t865 + t823 * t871;
t727 = t780 * t871 - t823 * t865;
t850 = legFrame(1,3);
t806 = sin(t850);
t809 = cos(t850);
t734 = t806 * t871 + t809 * t865;
t701 = t1026 * t734 + t724 * t809 + t727 * t806;
t874 = xDP(2);
t1004 = t701 * t874;
t731 = -t806 * t865 + t809 * t871;
t698 = t1026 * t731 - t724 * t806 + t727 * t809;
t875 = xDP(1);
t1007 = t698 * t875;
t812 = 0.1e1 / t823;
t1063 = 0.2e1 * pkin(1);
t864 = sin(qJ(2,1));
t750 = pkin(2) * t864 + pkin(3) * sin(t833);
t889 = 0.2e1 * qJ(2,1);
t832 = t889 + pkin(7);
t789 = sin(t832);
t972 = 0.2e1 * pkin(2) * pkin(3);
t894 = pkin(2) ^ 2;
t974 = t894 * sin(t889);
t777 = 0.2e1 * t833;
t892 = pkin(3) ^ 2;
t977 = t892 * sin(t777);
t907 = t789 * t972 + t974 + t977;
t710 = t1063 * t750 + t907;
t753 = t819 + t1026;
t744 = 0.1e1 / t753;
t992 = t744 * t873;
t947 = t710 * t992;
t677 = (t947 / 0.2e1 + t1004 + t1007) * t812;
t831 = qJ(2,2) + pkin(7);
t797 = cos(t831);
t1027 = pkin(3) * t797;
t818 = t868 * pkin(2);
t779 = t818 + pkin(1);
t822 = -pkin(6) + t852;
t863 = sin(qJ(1,2));
t869 = cos(qJ(1,2));
t723 = t779 * t863 + t822 * t869;
t726 = t779 * t869 - t822 * t863;
t849 = legFrame(2,3);
t805 = sin(t849);
t808 = cos(t849);
t733 = t805 * t869 + t808 * t863;
t700 = t1027 * t733 + t723 * t808 + t726 * t805;
t1005 = t700 * t874;
t730 = -t805 * t863 + t808 * t869;
t697 = t1027 * t730 - t723 * t805 + t726 * t808;
t1008 = t697 * t875;
t811 = 0.1e1 / t822;
t862 = sin(qJ(2,2));
t749 = pkin(2) * t862 + pkin(3) * sin(t831);
t887 = 0.2e1 * qJ(2,2);
t830 = t887 + pkin(7);
t787 = sin(t830);
t975 = t894 * sin(t887);
t776 = 0.2e1 * t831;
t978 = t892 * sin(t776);
t908 = t787 * t972 + t975 + t978;
t709 = t1063 * t749 + t908;
t752 = t818 + t1027;
t741 = 0.1e1 / t752;
t994 = t741 * t873;
t948 = t709 * t994;
t676 = (t948 / 0.2e1 + t1005 + t1008) * t811;
t829 = qJ(2,3) + pkin(7);
t794 = cos(t829);
t1028 = pkin(3) * t794;
t817 = t866 * pkin(2);
t778 = t817 + pkin(1);
t821 = -pkin(6) + t851;
t861 = sin(qJ(1,3));
t867 = cos(qJ(1,3));
t722 = t778 * t861 + t821 * t867;
t725 = t778 * t867 - t821 * t861;
t848 = legFrame(3,3);
t804 = sin(t848);
t807 = cos(t848);
t732 = t804 * t867 + t807 * t861;
t699 = t1028 * t732 + t722 * t807 + t725 * t804;
t1006 = t699 * t874;
t729 = -t804 * t861 + t807 * t867;
t696 = t1028 * t729 - t722 * t804 + t725 * t807;
t1009 = t696 * t875;
t810 = 0.1e1 / t821;
t860 = sin(qJ(2,3));
t748 = pkin(2) * t860 + pkin(3) * sin(t829);
t885 = 0.2e1 * qJ(2,3);
t828 = t885 + pkin(7);
t785 = sin(t828);
t976 = t894 * sin(t885);
t775 = 0.2e1 * t829;
t979 = t892 * sin(t775);
t909 = t785 * t972 + t976 + t979;
t708 = t1063 * t748 + t909;
t751 = t817 + t1028;
t738 = 0.1e1 / t751;
t996 = t738 * t873;
t949 = t708 * t996;
t675 = (t949 / 0.2e1 + t1006 + t1009) * t810;
t1062 = -0.2e1 * pkin(2);
t1061 = -t873 / 0.2e1;
t739 = 0.1e1 / t751 ^ 2;
t844 = t873 ^ 2;
t1060 = t739 * t844;
t742 = 0.1e1 / t752 ^ 2;
t1059 = t742 * t844;
t745 = 0.1e1 / t753 ^ 2;
t1058 = t745 * t844;
t1018 = Ifges(3,4) * t845;
t1020 = mrSges(3,2) * t845;
t1031 = m(3) * t894;
t1047 = 0.2e1 * mrSges(3,1);
t983 = t857 * t827;
t714 = 0.2e1 * t983 + (pkin(2) * t1047 + 0.4e1 * t1018) * t846 + t1031 + t1020 * t1062 + Ifges(2,2) - Ifges(2,1) - t857;
t1003 = t714 * t860;
t986 = t845 * t860;
t955 = pkin(3) * t986;
t1052 = -t955 + t1066;
t716 = 0.1e1 / t1052;
t1025 = pkin(3) * t845;
t719 = t1025 * t866 + t774 * t860;
t690 = (t716 * t719 * t873 + t729 * t875 + t732 * t874) * t810;
t841 = t866 ^ 2;
t1012 = t690 * t841;
t1017 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t1022 = (mrSges(2,3) + mrSges(3,3)) * pkin(5);
t1037 = mrSges(3,2) * pkin(1);
t1045 = -0.2e1 * t690;
t1046 = -2 * mrSges(3,3);
t1032 = -0.3e1 / 0.4e1 * t894;
t1035 = -t690 / 0.4e1;
t1029 = pkin(2) * t892;
t893 = pkin(2) * t894;
t1039 = -0.2e1 * t893 - 0.4e1 * t1029;
t784 = pkin(2) * t1024;
t968 = t784 + t894 / 0.2e1;
t1040 = -0.4e1 * pkin(1) * (t892 / 0.2e1 + t968);
t688 = pkin(1) * t690;
t667 = -t688 - (-t1009 / 0.2e1 - t1006 / 0.2e1 - t949 / 0.4e1) * t810;
t895 = pkin(1) ^ 2;
t681 = (t821 ^ 2 + t895) * t690;
t705 = t821 * t996;
t882 = 0.2e1 * pkin(7);
t791 = cos(t882 + qJ(2,3));
t795 = cos(qJ(2,3) - pkin(7));
t820 = t892 + t894;
t884 = 0.3e1 * qJ(2,3);
t891 = pkin(3) * t892;
t928 = 0.3e1 / 0.4e1 * t894;
t929 = 0.3e1 / 0.4e1 * t892;
t1023 = pkin(3) * t894;
t933 = -0.2e1 * t891 - 0.4e1 * t1023;
t997 = t738 * t810;
t939 = -t997 / 0.2e1;
t959 = -0.2e1 * t1023;
t960 = -0.2e1 * t1029;
t964 = -0.6e1 * t894 - 0.3e1 * t892;
t967 = cos(t828) + t846;
t657 = (t1039 * t866 + t791 * t960 + t794 * t933 + t795 * t959 + t1040) * t939 * t1060 + (t909 * t739 * t821 * t1061 + ((-t891 * cos(0.3e1 * t829) - t893 * cos(t884)) * t1035 - (t979 / 0.2e1 + t976 / 0.2e1) * t705 + (-(-cos(pkin(7) + t884) - t795) * t690 * t928 + (t1035 * t964 - t1063 * t675 + t681) * t794) * pkin(3) + ((-t1032 * t690 + t681) * t866 - (-cos(t882 + t884) - t791 - 0.2e1 * t866) * t690 * t929 + (-0.2e1 * t667 * t967 - t705 * t785) * pkin(3) - (pkin(3) * t967 + t1063 * t866) * t675) * pkin(2) + (-t675 / 0.2e1 - t667) * (cos(t775) * t892 + cos(t885) * t894 + t820)) * t738) * t810 * t690;
t1000 = t716 * t810;
t1041 = -0.2e1 * (t827 - 0.1e1 / 0.2e1) * t892 - 0.2e1 * t968;
t670 = -t688 + t675;
t982 = t892 * t690;
t991 = (0.2e1 * t784 + t820) * t844;
t660 = t739 * t810 / (t817 + (t846 * t866 - t986) * pkin(3)) * t991 - (t1012 * t1041 + t827 * t982 - t670 * t955 - t982 + t1052 * t675 - (t1045 * t955 - t670) * t1066) * t690 * t1000;
t880 = m(3) * pkin(2);
t813 = mrSges(2,1) + t880;
t1030 = pkin(1) * t813;
t684 = t690 * t1030;
t1016 = -t813 * pkin(5) + Ifges(2,5);
t934 = -m(3) * qJ(3,3) - mrSges(3,3);
t1036 = pkin(5) * mrSges(2,2);
t940 = Ifges(2,6) - t1036;
t693 = -(pkin(2) * t934 + t1016 + t1067) * t860 - t866 * (t1051 * t846 + t767 * t845 + t940);
t924 = mrSges(3,1) * t846 - t1020;
t736 = -t880 - t924;
t881 = m(3) * pkin(1);
t923 = mrSges(3,2) * t846 + t1021;
t711 = -t736 * t866 - t860 * t923 + t881;
t747 = 0.2e1 * t755;
t754 = (t813 - t1020) * t1063;
t879 = m(3) * pkin(5);
t770 = t879 - t934;
t883 = -pkin(5) / 0.2e1;
t814 = -qJ(3,3) / 0.2e1 + t883;
t840 = pkin(1) * t1047;
t877 = Ifges(3,6) / 0.2e1;
t878 = Ifges(3,5) / 0.2e1;
t905 = -m(2) * pkin(5) ^ 2 - (m(2) + m(3)) * t895 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + t983;
t906 = (mrSges(2,2) + t923) * t1063;
t943 = t738 * t748 * t1060;
t956 = 0.4e1 * t904;
t958 = 0.4e1 * t1019;
t961 = (mrSges(2,2) + t1021) * t1063;
t987 = (-t1036 / 0.2e1 + Ifges(2,6) / 0.2e1) * t873;
t971 = (-t714 * t841 - (t860 * t958 + (t747 * t860 + t840) * t846 + t754) * t866 + t860 * t961 - t851 ^ 2 * m(3) + qJ(3,3) * t1046 + t905) * t660 - t693 * t943 + t711 * t657 + 0.2e1 * (-t757 * t860 * t866 - (-t1037 * t860 - t1018) * t846 - t1022) * t660 - t738 * t956 * t1012 - ((pkin(2) * t770 + t1017 - t1067) * t996 - (t906 + 0.2e1 * t1003) * t690) * t866 * t996 - 0.2e1 * (((mrSges(3,2) * t814 + t877) * t996 - mrSges(3,1) * t688) * t846 + ((mrSges(3,1) * t814 + t878) * t996 + mrSges(3,2) * t688) * t845 + t738 * t987 - t684) * t860 * t996 + (-t675 * t770 - t738 * t904) * t1045;
t1002 = t714 * t862;
t985 = t845 * t862;
t954 = pkin(3) * t985;
t1053 = -t954 + t1065;
t717 = 0.1e1 / t1053;
t720 = t1025 * t868 + t774 * t862;
t691 = (t717 * t720 * t873 + t730 * t875 + t733 * t874) * t811;
t842 = t868 ^ 2;
t1011 = t691 * t842;
t1044 = -0.2e1 * t691;
t1034 = -t691 / 0.4e1;
t689 = pkin(1) * t691;
t668 = -t689 - (-t1008 / 0.2e1 - t1005 / 0.2e1 - t948 / 0.4e1) * t811;
t682 = (t822 ^ 2 + t895) * t691;
t706 = t822 * t994;
t792 = cos(t882 + qJ(2,2));
t798 = cos(qJ(2,2) - pkin(7));
t886 = 0.3e1 * qJ(2,2);
t995 = t741 * t811;
t938 = -t995 / 0.2e1;
t966 = cos(t830) + t846;
t658 = (t1039 * t868 + t792 * t960 + t797 * t933 + t798 * t959 + t1040) * t938 * t1059 + (t908 * t742 * t822 * t1061 + ((-t891 * cos(0.3e1 * t831) - t893 * cos(t886)) * t1034 - (t978 / 0.2e1 + t975 / 0.2e1) * t706 + (-(-cos(t886 + pkin(7)) - t798) * t691 * t928 + (t1034 * t964 - t1063 * t676 + t682) * t797) * pkin(3) + ((-t1032 * t691 + t682) * t868 - (-cos(t886 + t882) - t792 - 0.2e1 * t868) * t691 * t929 + (-0.2e1 * t668 * t966 - t706 * t787) * pkin(3) - (pkin(3) * t966 + t1063 * t868) * t676) * pkin(2) + (-t676 / 0.2e1 - t668) * (t892 * cos(t776) + t894 * cos(t887) + t820)) * t741) * t811 * t691;
t671 = -t689 + t676;
t981 = t892 * t691;
t999 = t717 * t811;
t661 = t742 * t811 / (t818 + (t846 * t868 - t985) * pkin(3)) * t991 - (t1011 * t1041 + t827 * t981 - t671 * t954 - t981 + t1053 * t676 - (t1044 * t954 - t671) * t1065) * t691 * t999;
t685 = t691 * t1030;
t935 = -m(3) * qJ(3,2) - mrSges(3,3);
t694 = -(pkin(2) * t935 + t1016 + t1068) * t862 - t868 * (t1050 * t846 + t768 * t845 + t940);
t712 = -t736 * t868 - t862 * t923 + t881;
t771 = t879 - t935;
t815 = -qJ(3,2) / 0.2e1 + t883;
t942 = t741 * t749 * t1059;
t970 = (-t714 * t842 - (t862 * t958 + (t747 * t862 + t840) * t846 + t754) * t868 + t862 * t961 - t852 ^ 2 * m(3) + qJ(3,2) * t1046 + t905) * t661 - t694 * t942 + t712 * t658 + 0.2e1 * (-t757 * t862 * t868 - (-t1037 * t862 - t1018) * t846 - t1022) * t661 - t741 * t956 * t1011 - ((pkin(2) * t771 + t1017 - t1068) * t994 - (t906 + 0.2e1 * t1002) * t691) * t868 * t994 - 0.2e1 * (((mrSges(3,2) * t815 + t877) * t994 - mrSges(3,1) * t689) * t846 + ((mrSges(3,1) * t815 + t878) * t994 + mrSges(3,2) * t689) * t845 + t741 * t987 - t685) * t862 * t994 + (-t676 * t771 - t741 * t904) * t1044;
t1001 = t714 * t864;
t984 = t845 * t864;
t953 = pkin(3) * t984;
t1054 = -t953 + t1064;
t718 = 0.1e1 / t1054;
t721 = t1025 * t870 + t774 * t864;
t692 = (t718 * t721 * t873 + t731 * t875 + t734 * t874) * t812;
t843 = t870 ^ 2;
t1010 = t692 * t843;
t1043 = -0.2e1 * t692;
t1033 = -t692 / 0.4e1;
t687 = t692 * pkin(1);
t666 = -t687 - (-t1007 / 0.2e1 - t1004 / 0.2e1 - t947 / 0.4e1) * t812;
t683 = (t823 ^ 2 + t895) * t692;
t707 = t823 * t992;
t800 = cos(qJ(2,1) + t882);
t802 = cos(qJ(2,1) - pkin(7));
t888 = 0.3e1 * qJ(2,1);
t993 = t744 * t812;
t937 = -t993 / 0.2e1;
t965 = cos(t832) + t846;
t659 = (t1039 * t870 + t800 * t960 + t801 * t933 + t802 * t959 + t1040) * t937 * t1058 + (t907 * t745 * t823 * t1061 + ((-t891 * cos(0.3e1 * t833) - t893 * cos(t888)) * t1033 - (t977 / 0.2e1 + t974 / 0.2e1) * t707 + (-(-cos(t888 + pkin(7)) - t802) * t692 * t928 + (t1033 * t964 - t1063 * t677 + t683) * t801) * pkin(3) + ((-t1032 * t692 + t683) * t870 - (-cos(t888 + t882) - t800 - 0.2e1 * t870) * t692 * t929 + (-0.2e1 * t666 * t965 - t707 * t789) * pkin(3) - (pkin(3) * t965 + t1063 * t870) * t677) * pkin(2) + (-t677 / 0.2e1 - t666) * (t892 * cos(t777) + t894 * cos(t889) + t820)) * t744) * t812 * t692;
t669 = -t687 + t677;
t980 = t892 * t692;
t998 = t718 * t812;
t662 = t745 * t812 / (t819 + (t846 * t870 - t984) * pkin(3)) * t991 - (t1010 * t1041 + t827 * t980 - t669 * t953 - t980 + t1054 * t677 - (t1043 * t953 - t669) * t1064) * t692 * t998;
t686 = t813 * t687;
t936 = -m(3) * qJ(3,1) - mrSges(3,3);
t695 = -(pkin(2) * t936 + t1016 + t1069) * t864 - t870 * (t1049 * t846 + t769 * t845 + t940);
t713 = -t736 * t870 - t864 * t923 + t881;
t772 = t879 - t936;
t816 = -qJ(3,1) / 0.2e1 + t883;
t941 = t744 * t750 * t1058;
t969 = (-t714 * t843 - (t864 * t958 + (t747 * t864 + t840) * t846 + t754) * t870 + t864 * t961 - t853 ^ 2 * m(3) + qJ(3,1) * t1046 + t905) * t662 - t695 * t941 + t713 * t659 + 0.2e1 * (-t757 * t864 * t870 - (-t1037 * t864 - t1018) * t846 - t1022) * t662 - t744 * t956 * t1010 - ((pkin(2) * t772 + t1017 - t1069) * t992 - (t906 + 0.2e1 * t1001) * t692) * t870 * t992 - 0.2e1 * (((mrSges(3,2) * t816 + t877) * t992 - mrSges(3,1) * t687) * t846 + ((mrSges(3,1) * t816 + t878) * t992 + mrSges(3,2) * t687) * t845 + t744 * t987 - t686) * t864 * t992 + (-t677 * t772 - t744 * t904) * t1043;
t1042 = -0.2e1 * t715;
t1038 = mrSges(2,2) * pkin(1);
t913 = t923 * t866;
t1015 = (t690 * (-m(3) * t851 + mrSges(3,3)) / 0.2e1 + (-t736 * t860 + t913) * t996) * t690;
t912 = t923 * t868;
t1014 = (t691 * (-m(3) * t852 + mrSges(3,3)) / 0.2e1 + (-t736 * t862 + t912) * t994) * t691;
t911 = t923 * t870;
t1013 = (t692 * (-m(3) * t853 + mrSges(3,3)) / 0.2e1 + (-t736 * t864 + t911) * t992) * t692;
t973 = -0.2e1 * t880;
t654 = -m(3) * t657 + t660 * t711;
t932 = t654 - 0.2e1 * t1015;
t655 = -m(3) * t658 + t661 * t712;
t931 = t655 - 0.2e1 * t1014;
t656 = -m(3) * t659 + t662 * t713;
t930 = t656 - 0.2e1 * t1013;
t728 = t924 * t1062 - Ifges(2,3) - Ifges(3,3) - t1031;
t1 = [-(t698 * t930 + t731 * t969) * t812 - (t697 * t931 + t730 * t970) * t811 - (t696 * t932 + t729 * t971) * t810; -(t701 * t930 + t734 * t969) * t812 - (t700 * t931 + t733 * t970) * t811 - (t699 * t932 + t732 * t971) * t810; -t971 * t719 * t1000 - t970 * t720 * t999 - t969 * t721 * t998 + (t662 * t695 - t728 * t941 - ((-t677 * t973 - t686) * t864 + (t864 * t924 + t911) * (-t687 - (-0.2e1 * t1004 - t947 - 0.2e1 * t1007) * t812) - (t843 * t1042 + (t1001 + t1038) * t870 + t715) * t692) * t692) * t744 + (t661 * t694 - t728 * t942 - ((-t676 * t973 - t685) * t862 + (t862 * t924 + t912) * (-t689 - (-0.2e1 * t1005 - t948 - 0.2e1 * t1008) * t811) - (t842 * t1042 + (t1002 + t1038) * t868 + t715) * t691) * t691) * t741 + (t660 * t693 - t728 * t943 - ((-t675 * t973 - t684) * t860 + (t860 * t924 + t913) * (-t688 - (-0.2e1 * t1006 - t949 - 0.2e1 * t1009) * t810) - (t841 * t1042 + (t1003 + t1038) * t866 + t715) * t690) * t690) * t738 + (t1013 * t993 + t656 * t937) * t710 + (t1014 * t995 + t655 * t938) * t709 + (t1015 * t997 + t654 * t939) * t708;];
taucX  = t1;
