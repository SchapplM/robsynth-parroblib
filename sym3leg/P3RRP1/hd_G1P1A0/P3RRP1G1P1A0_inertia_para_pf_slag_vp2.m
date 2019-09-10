% Calculate inertia matrix for parallel robot
% P3RRP1G1P1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:31
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRP1G1P1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1G1P1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:05:28
% EndTime: 2019-05-03 15:05:31
% DurationCPUTime: 3.67s
% Computational Cost: add. (23150->452), mult. (50282->640), div. (432->3), fcn. (21104->26), ass. (0->286)
t897 = pkin(1) ^ 2;
t987 = Ifges(1,3) + (m(2) + m(3)) * t897;
t887 = (qJ(3,3) ^ 2);
t879 = cos(qJ(1,3));
t859 = t879 ^ 2;
t935 = t859 * t897;
t986 = (2 * t887) + t897 - t935;
t888 = (qJ(3,2) ^ 2);
t881 = cos(qJ(1,2));
t860 = t881 ^ 2;
t934 = t860 * t897;
t985 = (2 * t888) + t897 - t934;
t889 = (qJ(3,1) ^ 2);
t883 = cos(qJ(1,1));
t861 = t883 ^ 2;
t933 = t861 * t897;
t984 = (2 * t889) + t897 - t933;
t983 = 2 * pkin(1);
t982 = 2 * mrSges(3,3);
t869 = legFrame(3,3);
t836 = sin(t869);
t811 = t836 * qJ(3,3);
t807 = pkin(2) * t811;
t981 = -0.2e1 * t807;
t870 = legFrame(2,3);
t837 = sin(t870);
t812 = t837 * qJ(3,2);
t808 = pkin(2) * t812;
t980 = -0.2e1 * t808;
t871 = legFrame(1,3);
t838 = sin(t871);
t813 = t838 * qJ(3,1);
t809 = pkin(2) * t813;
t979 = -0.2e1 * t809;
t873 = sin(qJ(1,3));
t978 = -0.2e1 * t873;
t875 = sin(qJ(1,2));
t977 = -0.2e1 * t875;
t877 = sin(qJ(1,1));
t976 = -0.2e1 * t877;
t975 = -0.2e1 * t879;
t974 = -0.2e1 * t881;
t973 = -0.2e1 * t883;
t864 = qJ(1,1) + qJ(2,1);
t832 = sin(t864);
t972 = pkin(1) * t832;
t971 = pkin(2) * qJ(3,1);
t970 = pkin(2) * qJ(3,2);
t969 = pkin(2) * qJ(3,3);
t839 = cos(t869);
t968 = pkin(2) * t839;
t840 = cos(t870);
t967 = pkin(2) * t840;
t841 = cos(t871);
t966 = pkin(2) * t841;
t965 = pkin(2) * t879;
t964 = pkin(2) * t881;
t963 = pkin(2) * t883;
t814 = t836 * pkin(2);
t815 = t837 * pkin(2);
t816 = t838 * pkin(2);
t962 = mrSges(3,3) - mrSges(2,2);
t961 = qJ(3,1) * t841;
t960 = qJ(3,1) * t877;
t959 = qJ(3,2) * t840;
t958 = qJ(3,2) * t875;
t957 = qJ(3,3) * t839;
t956 = qJ(3,3) * t873;
t780 = t811 + t968;
t955 = t780 * t873;
t784 = t812 + t967;
t954 = t784 * t875;
t788 = t813 + t966;
t953 = t788 * t877;
t862 = qJ(1,3) + qJ(2,3);
t833 = cos(t862);
t952 = t811 * t833;
t863 = qJ(1,2) + qJ(2,2);
t834 = cos(t863);
t951 = t812 * t834;
t835 = cos(t864);
t950 = t813 * t835;
t896 = pkin(2) ^ 2;
t843 = -t887 + t896;
t827 = 0.1e1 + t843;
t949 = t827 * t839;
t844 = -t888 + t896;
t828 = 0.1e1 + t844;
t948 = t828 * t840;
t845 = -t889 + t896;
t829 = 0.1e1 + t845;
t947 = t829 * t841;
t830 = sin(t862);
t946 = t830 * qJ(3,3);
t945 = t830 * t833;
t944 = t830 * t836;
t831 = sin(t863);
t943 = t831 * qJ(3,2);
t942 = t831 * t834;
t941 = t831 * t837;
t940 = t832 * t835;
t939 = t832 * t838;
t938 = t833 * t839;
t937 = t834 * t840;
t936 = t835 * t841;
t932 = t873 * t879;
t931 = t875 * t881;
t930 = t877 * t883;
t929 = -0.1e1 / 0.2e1 - t896 / 0.2e1;
t867 = t896 / 0.2e1;
t928 = 0.1e1 / 0.2e1 + t867;
t927 = pkin(2) * t961;
t926 = pkin(2) * t959;
t925 = pkin(2) * t957;
t924 = 0.2e1 * pkin(2) * mrSges(3,1) + Ifges(3,2) + Ifges(2,3);
t842 = -m(3) * pkin(2) - mrSges(3,1);
t923 = t839 * t946;
t922 = t840 * t943;
t921 = t832 * t961;
t920 = t897 * t932;
t919 = t897 * t931;
t918 = t897 * t930;
t884 = xP(3);
t850 = sin(t884);
t851 = cos(t884);
t885 = mrSges(4,2);
t886 = mrSges(4,1);
t917 = -t850 * t885 + t851 * t886;
t916 = -(2 * t887) - t935;
t915 = -(2 * t888) - t934;
t914 = -(2 * t889) - t933;
t913 = t839 * t920;
t912 = t840 * t919;
t911 = t841 * t918;
t771 = (t887 + t896) * m(3) + qJ(3,3) * t982 + t924;
t772 = (t888 + t896) * m(3) + qJ(3,2) * t982 + t924;
t773 = (t889 + t896) * m(3) + qJ(3,1) * t982 + t924;
t910 = -t850 * t886 - t851 * t885;
t781 = t811 - t968;
t782 = t814 + t957;
t909 = -t781 * t859 - t782 * t932;
t908 = -t781 * t932 + t782 * t859;
t785 = t812 - t967;
t786 = t815 + t959;
t907 = -t785 * t860 - t786 * t931;
t906 = -t785 * t931 + t786 * t860;
t789 = t813 - t966;
t790 = t816 + t961;
t905 = -t789 * t861 - t790 * t930;
t904 = -t789 * t930 + t790 * t861;
t823 = mrSges(2,1) - t842;
t878 = cos(qJ(2,3));
t903 = ((m(3) * qJ(3,3) + t962) * sin(qJ(2,3)) + t823 * t878) * pkin(1);
t880 = cos(qJ(2,2));
t902 = ((m(3) * qJ(3,2) + t962) * sin(qJ(2,2)) + t823 * t880) * pkin(1);
t882 = cos(qJ(2,1));
t901 = ((m(3) * qJ(3,1) + t962) * sin(qJ(2,1)) + t823 * t882) * pkin(1);
t900 = -t829 * t861 + 0.2e1 * t930 * t971;
t899 = -t828 * t860 + 0.2e1 * t931 * t970;
t898 = -t827 * t859 + 0.2e1 * t932 * t969;
t895 = koppelP(1,1);
t894 = koppelP(2,1);
t893 = koppelP(3,1);
t892 = koppelP(1,2);
t891 = koppelP(2,2);
t890 = koppelP(3,2);
t866 = 0.1e1 + t897;
t865 = 0.2e1 + t897;
t858 = -t889 / 0.2e1;
t857 = -t888 / 0.2e1;
t856 = -t887 / 0.2e1;
t826 = t835 ^ 2;
t825 = t834 ^ 2;
t824 = t833 ^ 2;
t806 = -0.2e1 * t927;
t805 = -0.2e1 * t926;
t804 = -0.2e1 * t925;
t803 = 0.2e1 * t809;
t802 = 0.2e1 * t808;
t801 = 0.2e1 * t807;
t800 = -mrSges(3,1) + (-pkin(1) * t882 - pkin(2)) * m(3);
t799 = -mrSges(3,1) + (-pkin(1) * t880 - pkin(2)) * m(3);
t798 = -mrSges(3,1) + (-pkin(1) * t878 - pkin(2)) * m(3);
t797 = -t960 + t963;
t796 = -t958 + t964;
t795 = -t956 + t965;
t794 = t838 * t918;
t793 = t837 * t919;
t792 = t836 * t920;
t791 = t816 - t961;
t787 = t815 - t959;
t783 = t814 - t957;
t779 = -t850 * t892 + t851 * t895;
t778 = -t850 * t891 + t851 * t894;
t777 = -t850 * t890 + t851 * t893;
t776 = -t850 * t895 - t851 * t892;
t775 = -t850 * t894 - t851 * t891;
t774 = -t850 * t893 - t851 * t890;
t770 = t838 * t845 + t806;
t769 = t837 * t844 + t805;
t768 = t836 * t843 + t804;
t767 = t829 * t838 + t806;
t766 = t828 * t837 + t805;
t765 = t827 * t836 + t804;
t764 = t767 * t883;
t763 = t766 * t881;
t762 = t765 * t879;
t761 = t791 * t883 - t953;
t760 = t788 * t883 + t791 * t877;
t759 = t787 * t881 - t954;
t758 = t784 * t881 + t787 * t875;
t757 = t783 * t879 - t955;
t756 = t780 * t879 + t783 * t873;
t755 = t901 + t773;
t754 = t902 + t772;
t753 = t903 + t771;
t752 = t773 + 0.2e1 * t901 + t987;
t751 = t772 + 0.2e1 * t902 + t987;
t750 = t771 + 0.2e1 * t903 + t987;
t749 = (t841 * t845 + t803) * t883 + t770 * t877;
t748 = (t840 * t844 + t802) * t881 + t769 * t875;
t747 = (t839 * t843 + t801) * t879 + t768 * t873;
t746 = (t803 + t947) * t883 + t877 * t767;
t745 = (t802 + t948) * t881 + t875 * t766;
t744 = (t801 + t949) * t879 + t873 * t765;
t743 = t770 * t883 + ((t867 + t858) * t841 + t809) * t976;
t742 = t769 * t881 + ((t867 + t857) * t840 + t808) * t977;
t741 = t768 * t879 + ((t867 + t856) * t839 + t807) * t978;
t740 = 0.1e1 / (((t829 * t930 + (0.2e1 * t861 - 0.1e1) * t971) * t972 + t960) * t835 * t983 + qJ(3,1) * t972 * t973 - t889 * t865 + (0.2e1 * (t889 / 0.2e1 - t900 + t929) * t826 + t900) * t897);
t739 = 0.1e1 / ((pkin(1) * (t828 * t931 + (0.2e1 * t860 - 0.1e1) * t970) * t831 + t958) * t834 * t983 + pkin(1) * t943 * t974 - t888 * t865 + (0.2e1 * (t888 / 0.2e1 - t899 + t929) * t825 + t899) * t897);
t738 = 0.1e1 / ((pkin(1) * (t827 * t932 + (0.2e1 * t859 - 0.1e1) * t969) * t830 + t956) * t833 * t983 + pkin(1) * t946 * t975 - t887 * t865 + (0.2e1 * (t887 / 0.2e1 - t898 + t929) * t824 + t898) * t897);
t737 = (-t936 + t939) * qJ(3,1) + (-t764 * t826 - t746 * t940 + (t838 * t896 + t838 - t927) * t883 + (-(t979 - t947) * t826 - qJ(3,1) * t791) * t877) * pkin(1);
t736 = (-t937 + t941) * qJ(3,2) + (-t763 * t825 - t745 * t942 + (t837 * t896 + t837 - t926) * t881 + (-(t980 - t948) * t825 - qJ(3,2) * t787) * t875) * pkin(1);
t735 = (-t938 + t944) * qJ(3,3) + (-t762 * t824 - t744 * t945 + (t836 * t896 + t836 - t925) * t879 + (-(t981 - t949) * t824 - qJ(3,3) * t783) * t873) * pkin(1);
t734 = -t921 - t950 + (t746 * t826 - (t764 + ((t858 + t928) * t841 + t809) * t976) * t940 - (t841 * t896 + t809 + t841) * t883 + qJ(3,1) * t953) * pkin(1);
t733 = -t922 - t951 + (t745 * t825 - (t763 + ((t857 + t928) * t840 + t808) * t977) * t942 - (t840 * t896 + t808 + t840) * t881 + qJ(3,2) * t954) * pkin(1);
t732 = -t923 - t952 + (t744 * t824 - (t762 + ((t856 + t928) * t839 + t807) * t978) * t945 - (t839 * t896 + t807 + t839) * t879 + qJ(3,3) * t955) * pkin(1);
t731 = (-t838 * t984 + t806 + t911) * t835 + (t914 * t841 + t794 + t803) * t832 + (-t761 * t826 - t760 * t940 + (t816 - 0.2e1 * t961) * t883 + t877 * t813) * pkin(1);
t730 = (-t837 * t985 + t805 + t912) * t834 + (t915 * t840 + t793 + t802) * t831 + (-t759 * t825 - t758 * t942 + (t815 - 0.2e1 * t959) * t881 + t875 * t812) * pkin(1);
t729 = (-t836 * t986 + t804 + t913) * t833 + (t916 * t839 + t792 + t801) * t830 + (-t757 * t824 - t756 * t945 + (t814 - 0.2e1 * t957) * t879 + t873 * t811) * pkin(1);
t728 = (t841 * t984 + t794 + t979) * t835 + (t914 * t838 + t806 - t911) * t832 + (-t761 * t940 + t760 * t826 + t813 * t973 + (-t960 - t963) * t841) * pkin(1);
t727 = (t840 * t985 + t793 + t980) * t834 + (t915 * t837 + t805 - t912) * t831 + (-t759 * t942 + t758 * t825 + t812 * t974 + (-t958 - t964) * t840) * pkin(1);
t726 = (t839 * t986 + t792 + t981) * t833 + (t916 * t836 + t804 - t913) * t830 + (-t757 * t945 + t756 * t824 + t811 * t975 + (-t956 - t965) * t839) * pkin(1);
t725 = -t866 * t921 - t950 + (t743 * t940 - t749 * t826 + t788 * t797) * pkin(1) + ((t905 - t966) * t835 + t904 * t832) * t897;
t724 = -t866 * t922 - t951 + (t742 * t942 - t748 * t825 + t784 * t796) * pkin(1) + ((t907 - t967) * t834 + t906 * t831) * t897;
t723 = -t866 * t923 - t952 + (t741 * t945 - t747 * t824 + t780 * t795) * pkin(1) + ((t909 - t968) * t833 + t908 * t830) * t897;
t722 = (t866 * t939 - t936) * qJ(3,1) + (t743 * t826 + t749 * t940 - t791 * t797) * pkin(1) + ((-t904 + t816) * t835 + t905 * t832) * t897;
t721 = (t866 * t941 - t937) * qJ(3,2) + (t742 * t825 + t748 * t942 - t787 * t796) * pkin(1) + ((-t906 + t815) * t834 + t907 * t831) * t897;
t720 = (t866 * t944 - t938) * qJ(3,3) + (t741 * t824 + t747 * t945 - t783 * t795) * pkin(1) + ((-t908 + t814) * t833 + t909 * t830) * t897;
t719 = (t734 * t779 + t737 * t776) * t740;
t718 = (t733 * t778 + t736 * t775) * t739;
t717 = (t732 * t777 + t735 * t774) * t738;
t716 = (t728 * t779 + t731 * t776) * t740;
t715 = (t727 * t778 + t730 * t775) * t739;
t714 = (t726 * t777 + t729 * t774) * t738;
t713 = (t722 * t776 + t725 * t779) * t740;
t712 = (t721 * t775 + t724 * t778) * t739;
t711 = (t720 * t774 + t723 * t777) * t738;
t710 = (m(3) * t731 + t722 * t842 + t737 * t800) * t740;
t709 = (m(3) * t730 + t721 * t842 + t736 * t799) * t739;
t708 = (m(3) * t729 + t720 * t842 + t735 * t798) * t738;
t707 = (m(3) * t728 + t725 * t842 + t734 * t800) * t740;
t706 = (m(3) * t727 + t724 * t842 + t733 * t799) * t739;
t705 = (m(3) * t726 + t723 * t842 + t732 * t798) * t738;
t704 = (t722 * t773 + t731 * t842 + t737 * t755) * t740;
t703 = (t721 * t772 + t730 * t842 + t736 * t754) * t739;
t702 = (t720 * t771 + t729 * t842 + t735 * t753) * t738;
t701 = (t725 * t773 + t728 * t842 + t734 * t755) * t740;
t700 = (t724 * t772 + t727 * t842 + t733 * t754) * t739;
t699 = (t723 * t771 + t726 * t842 + t732 * t753) * t738;
t698 = (t722 * t755 + t731 * t800 + t737 * t752) * t740;
t697 = (t721 * t754 + t730 * t799 + t736 * t751) * t739;
t696 = (t720 * t753 + t729 * t798 + t735 * t750) * t738;
t695 = (t725 * t755 + t728 * t800 + t734 * t752) * t740;
t694 = (t724 * t754 + t727 * t799 + t733 * t751) * t739;
t693 = (t723 * t753 + t726 * t798 + t732 * t750) * t738;
t692 = t716 * m(3) + t713 * t842 + t719 * t800;
t691 = t715 * m(3) + t712 * t842 + t718 * t799;
t690 = t714 * m(3) + t711 * t842 + t717 * t798;
t689 = t713 * t773 + t716 * t842 + t719 * t755;
t688 = t712 * t772 + t715 * t842 + t718 * t754;
t687 = t711 * t771 + t714 * t842 + t717 * t753;
t686 = t713 * t755 + t716 * t800 + t719 * t752;
t685 = t712 * t754 + t715 * t799 + t718 * t751;
t684 = t711 * t753 + t714 * t798 + t717 * t750;
t1 = [m(4) + (t698 * t737 + t704 * t722 + t710 * t731) * t740 + (t697 * t736 + t703 * t721 + t709 * t730) * t739 + (t696 * t735 + t702 * t720 + t708 * t729) * t738, (t698 * t734 + t704 * t725 + t710 * t728) * t740 + (t697 * t733 + t703 * t724 + t709 * t727) * t739 + (t696 * t732 + t702 * t723 + t708 * t726) * t738, t696 * t717 + t697 * t718 + t698 * t719 + t702 * t711 + t703 * t712 + t704 * t713 + t708 * t714 + t709 * t715 + t710 * t716 + t910; (t695 * t737 + t701 * t722 + t707 * t731) * t740 + (t694 * t736 + t700 * t721 + t706 * t730) * t739 + (t693 * t735 + t699 * t720 + t705 * t729) * t738, m(4) + (t695 * t734 + t701 * t725 + t707 * t728) * t740 + (t694 * t733 + t700 * t724 + t706 * t727) * t739 + (t693 * t732 + t699 * t723 + t705 * t726) * t738, t693 * t717 + t694 * t718 + t695 * t719 + t699 * t711 + t700 * t712 + t701 * t713 + t705 * t714 + t706 * t715 + t707 * t716 + t917; (t686 * t737 + t689 * t722 + t692 * t731) * t740 + (t685 * t736 + t688 * t721 + t691 * t730) * t739 + (t684 * t735 + t687 * t720 + t690 * t729) * t738 + t910, (t686 * t734 + t689 * t725 + t692 * t728) * t740 + (t685 * t733 + t688 * t724 + t691 * t727) * t739 + (t684 * t732 + t687 * t723 + t690 * t726) * t738 + t917, t684 * t717 + t685 * t718 + t686 * t719 + t687 * t711 + t688 * t712 + t689 * t713 + t690 * t714 + t691 * t715 + t692 * t716 + Ifges(4,3);];
MX  = t1;
