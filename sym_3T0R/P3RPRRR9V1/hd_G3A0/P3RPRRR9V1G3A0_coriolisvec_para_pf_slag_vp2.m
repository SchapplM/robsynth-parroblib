% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR9V1G3A0
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:37
% EndTime: 2020-08-06 18:56:42
% DurationCPUTime: 5.10s
% Computational Cost: add. (20541->458), mult. (26184->747), div. (4311->10), fcn. (22560->36), ass. (0->304)
t1016 = -2 * pkin(5);
t831 = sin(pkin(7));
t1010 = -0.2e1 * t831;
t832 = cos(pkin(7));
t790 = t832 * pkin(2);
t1015 = -0.2e1 * t790;
t1013 = 2 * pkin(1);
t1014 = t1013 / 0.2e1;
t822 = pkin(7) + qJ(3,1);
t789 = cos(t822);
t1002 = pkin(3) * t789;
t777 = t790 + pkin(1);
t752 = t777 + t1002;
t821 = pkin(7) + qJ(3,2);
t788 = cos(t821);
t1003 = pkin(3) * t788;
t751 = t777 + t1003;
t820 = pkin(7) + qJ(3,3);
t787 = cos(t820);
t1004 = pkin(3) * t787;
t750 = t777 + t1004;
t833 = pkin(5) + qJ(2,3);
t768 = t833 * mrSges(3,1) - Ifges(3,5);
t834 = pkin(5) + qJ(2,2);
t769 = t834 * mrSges(3,1) - Ifges(3,5);
t835 = pkin(5) + qJ(2,1);
t770 = t835 * mrSges(3,1) - Ifges(3,5);
t765 = t833 * mrSges(3,2) - Ifges(3,6);
t766 = t834 * mrSges(3,2) - Ifges(3,6);
t767 = t835 * mrSges(3,2) - Ifges(3,6);
t1012 = 0.2e1 * pkin(3);
t819 = t832 ^ 2;
t1011 = 0.2e1 * t819;
t1009 = 0.4e1 * t832;
t1008 = -4 * pkin(5) - 4 * pkin(6);
t1007 = pkin(1) * mrSges(3,2);
t1006 = pkin(2) * mrSges(3,1);
t1005 = pkin(2) * mrSges(3,2);
t784 = sin(t820);
t838 = legFrame(3,2);
t804 = cos(t838);
t801 = sin(t838);
t848 = cos(qJ(1,3));
t967 = t801 * t848;
t738 = t804 * t784 - t787 * t967;
t964 = t804 * t848;
t739 = t801 * t784 + t787 * t964;
t781 = 0.1e1 / t787;
t816 = pkin(6) + t833;
t798 = 1 / t816;
t842 = sin(qJ(1,3));
t854 = xDP(3);
t855 = xDP(2);
t856 = xDP(1);
t711 = (-t842 * t854 + (t738 * t855 + t739 * t856) * t781) * t798;
t709 = pkin(1) * t711;
t785 = sin(t821);
t839 = legFrame(2,2);
t805 = cos(t839);
t802 = sin(t839);
t850 = cos(qJ(1,2));
t966 = t802 * t850;
t740 = t805 * t785 - t788 * t966;
t963 = t805 * t850;
t741 = t802 * t785 + t788 * t963;
t782 = 0.1e1 / t788;
t817 = pkin(6) + t834;
t799 = 1 / t817;
t844 = sin(qJ(1,2));
t712 = (-t844 * t854 + (t740 * t855 + t741 * t856) * t782) * t799;
t710 = pkin(1) * t712;
t791 = pkin(1) * t831;
t786 = sin(t822);
t840 = legFrame(1,2);
t806 = cos(t840);
t803 = sin(t840);
t852 = cos(qJ(1,1));
t965 = t803 * t852;
t742 = t806 * t786 - t789 * t965;
t962 = t806 * t852;
t743 = t803 * t786 + t789 * t962;
t783 = 0.1e1 / t789;
t818 = pkin(6) + t835;
t800 = 1 / t818;
t846 = sin(qJ(1,1));
t713 = (-t846 * t854 + (t742 * t855 + t743 * t856) * t783) * t800;
t708 = t713 * pkin(1);
t847 = cos(qJ(3,3));
t827 = t847 ^ 2;
t1001 = t827 * pkin(3);
t849 = cos(qJ(3,2));
t828 = t849 ^ 2;
t1000 = t828 * pkin(3);
t851 = cos(qJ(3,1));
t829 = t851 ^ 2;
t999 = t829 * pkin(3);
t998 = t847 * pkin(2);
t997 = t849 * pkin(2);
t996 = t851 * pkin(2);
t837 = mrSges(2,3) + mrSges(3,3);
t836 = Ifges(3,1) - Ifges(3,2);
t995 = Ifges(2,4) - Ifges(3,4);
t994 = mrSges(3,1) * t831;
t993 = mrSges(3,2) * t831;
t841 = sin(qJ(3,3));
t992 = mrSges(3,2) * t841;
t843 = sin(qJ(3,2));
t991 = mrSges(3,2) * t843;
t845 = sin(qJ(3,1));
t990 = mrSges(3,2) * t845;
t795 = mrSges(3,2) * t847;
t796 = mrSges(3,2) * t849;
t797 = mrSges(3,2) * t851;
t989 = Ifges(3,4) * t841;
t988 = Ifges(3,4) * t843;
t987 = Ifges(3,4) * t845;
t986 = t827 * Ifges(3,4);
t985 = t828 * Ifges(3,4);
t984 = t829 * Ifges(3,4);
t792 = t841 * mrSges(3,1);
t793 = t843 * mrSges(3,1);
t794 = t845 * mrSges(3,1);
t983 = m(3) * pkin(2) + mrSges(2,1);
t982 = t711 * t819;
t981 = t711 * t831;
t980 = t712 * t819;
t979 = t712 * t831;
t978 = t713 * t819;
t977 = t713 * t831;
t877 = 0.1e1 / pkin(3);
t729 = (t801 * t856 + t804 * t855) * t877 * t781;
t976 = t729 ^ 2 * t781;
t730 = (t802 * t856 + t805 * t855) * t877 * t782;
t975 = t730 ^ 2 * t782;
t731 = (t803 * t856 + t806 * t855) * t877 * t783;
t974 = t731 ^ 2 * t783;
t961 = t831 * t841;
t747 = 0.1e1 / (t832 * t847 - t961);
t973 = t747 * t798;
t960 = t831 * t843;
t748 = 0.1e1 / (t832 * t849 - t960);
t972 = t748 * t799;
t959 = t831 * t845;
t749 = 0.1e1 / (t832 * t851 - t959);
t971 = t749 * t800;
t872 = pkin(2) / 0.2e1;
t970 = (t847 * pkin(3) + t872) * t841;
t969 = (t849 * pkin(3) + t872) * t843;
t968 = (t851 * pkin(3) + t872) * t845;
t958 = t836 * t827;
t957 = t836 * t828;
t956 = t836 * t829;
t955 = t836 * t841;
t954 = t836 * t843;
t953 = t836 * t845;
t937 = pkin(1) * t848 + t842 * t816;
t723 = t937 * t961 + (t827 - 0.1e1) * t848 * pkin(3);
t735 = pkin(1) * t841 + (-pkin(3) + t998 + 0.2e1 * t1001) * t831;
t871 = -pkin(3) / 0.2e1;
t753 = t1001 + t998 / 0.2e1 + t871;
t907 = t848 * t961;
t883 = pkin(2) * t907 + (t907 * t1012 - t937) * t847;
t946 = t847 * (-t841 * pkin(3) + t791);
t699 = (-t753 * t967 + t804 * t970) * t1011 + (t804 * t735 + t883 * t801) * t832 + t723 * t801 + t804 * t946;
t690 = t699 * t855 * t973;
t700 = (t753 * t964 + t801 * t970) * t1011 + (t801 * t735 - t883 * t804) * t832 - t723 * t804 + t801 * t946;
t691 = t700 * t856 * t973;
t732 = -t750 * t842 + t816 * t848;
t720 = t732 * t798 * t854;
t685 = t709 - 0.2e1 * t691 - 0.2e1 * t690 - 0.2e1 * t720;
t952 = t841 * t685;
t951 = t841 * t847;
t936 = pkin(1) * t850 + t844 * t817;
t724 = t936 * t960 + (t828 - 0.1e1) * t850 * pkin(3);
t736 = pkin(1) * t843 + (-pkin(3) + t997 + 0.2e1 * t1000) * t831;
t754 = t1000 + t997 / 0.2e1 + t871;
t906 = t850 * t960;
t882 = pkin(2) * t906 + (t906 * t1012 - t936) * t849;
t945 = t849 * (-t843 * pkin(3) + t791);
t701 = (-t754 * t966 + t805 * t969) * t1011 + (t805 * t736 + t882 * t802) * t832 + t724 * t802 + t805 * t945;
t692 = t701 * t855 * t972;
t702 = (t754 * t963 + t802 * t969) * t1011 + (t802 * t736 - t882 * t805) * t832 - t724 * t805 + t802 * t945;
t693 = t702 * t856 * t972;
t733 = -t751 * t844 + t817 * t850;
t721 = t733 * t799 * t854;
t686 = t710 - 0.2e1 * t693 - 0.2e1 * t692 - 0.2e1 * t721;
t950 = t843 * t686;
t949 = t843 * t849;
t935 = pkin(1) * t852 + t846 * t818;
t725 = t935 * t959 + (t829 - 0.1e1) * t852 * pkin(3);
t737 = pkin(1) * t845 + (-pkin(3) + t996 + 0.2e1 * t999) * t831;
t755 = t999 + t996 / 0.2e1 + t871;
t905 = t852 * t959;
t881 = pkin(2) * t905 + (t905 * t1012 - t935) * t851;
t944 = t851 * (-t845 * pkin(3) + t791);
t703 = (-t755 * t965 + t806 * t968) * t1011 + (t806 * t737 + t881 * t803) * t832 + t725 * t803 + t806 * t944;
t694 = t703 * t855 * t971;
t704 = (t755 * t962 + t803 * t968) * t1011 + (t803 * t737 - t881 * t806) * t832 - t725 * t806 + t803 * t944;
t695 = t704 * t856 * t971;
t734 = -t752 * t846 + t818 * t852;
t722 = t734 * t800 * t854;
t684 = t708 - 0.2e1 * t695 - 0.2e1 * t694 - 0.2e1 * t722;
t948 = t845 * t684;
t947 = t845 * t851;
t682 = t709 - t691 / 0.2e1 - t690 / 0.2e1 - t720 / 0.2e1;
t687 = t691 + t690 + t720;
t778 = 0.2e1 * t820;
t868 = 0.2e1 * pkin(7);
t873 = qJ(2,3) ^ 2;
t876 = pkin(3) ^ 2;
t878 = pkin(2) ^ 2;
t879 = pkin(1) ^ 2;
t880 = -t878 * cos(t868) - (2 * pkin(6) ^ 2) - t876 - t878 - (2 * t879) + ((-4 * pkin(6) + t1016) * pkin(5));
t669 = ((t682 * t1015 + ((qJ(2,3) * t1008) - t876 * cos(t778) - (2 * t873) + t880) * t711 / 0.2e1 + (t750 + t1014) * t687) * t711 + ((t729 * t816 * t784 - 0.2e1 * t682 * t787 + (-cos(qJ(3,3) + t868) * pkin(2) - t998) * t711) * t711 - (-t711 * t816 * sin(t778) / 0.2e1 + t729 * t750) * t781 * t729) * pkin(3)) * t798;
t678 = (-pkin(3) * t976 + (-t709 + 0.2e1 * t687 + (-t790 - t1004) * t711) * t711) * t798;
t705 = t711 * t955;
t714 = -(-t765 * t847 - t768 * t841) * t832 - t831 * (t841 * t765 - t768 * t847);
t858 = m(2) + m(3);
t814 = pkin(1) * t858;
t893 = mrSges(3,1) * t847 - t992;
t940 = t795 + t792;
t717 = -(-t893 - t983) * t832 - (mrSges(2,2) + t940) * t831 + t814;
t917 = m(3) * pkin(5) + t837;
t756 = t858 * qJ(2,3) + t917;
t812 = -t1005 / 0.2e1;
t813 = t1006 / 0.4e1;
t824 = 0.2e1 * t1006;
t859 = Ifges(3,6) / 0.2e1;
t860 = Ifges(3,5) / 0.4e1;
t869 = -pkin(5) / 0.4e1;
t870 = -pkin(5) / 0.2e1;
t884 = (mrSges(3,3) * t1016) - t858 * t879 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3);
t894 = t878 * m(3) - Ifges(2,1) + Ifges(2,2) + t836;
t897 = t958 * t981;
t898 = (mrSges(3,1) * t1013) + t1005 * t1010;
t910 = t784 * t976;
t923 = pkin(2) * t992;
t913 = (t923 - t836) * t981;
t916 = t711 * t986;
t919 = -t1006 / 0.2e1;
t920 = -t1007 / 0.2e1;
t926 = mrSges(3,1) * t709;
t930 = 0.4e1 * t989;
t931 = 0.4e1 * Ifges(3,4) * t831;
t933 = 0.2e1 * t791;
t934 = -2 * t1007;
t943 = (-((t824 + t930) * t847 - 0.2e1 * t923 + t894) * t819 - (t827 * t931 + t898 * t847 + (t983 - t992) * t1013) * t832 - t958 + (t792 + mrSges(2,2)) * t933 - (t833 ^ 2) * m(3) - m(2) * t873 + t884) * t678 + t717 * t669 - t714 * t910 + 0.2e1 * (Ifges(3,4) * t951 + t819 * t958 - t837 * qJ(2,3) + (-(-pkin(2) * t792 + t836 * t951 + t995) * t832 + pkin(1) * t795) * t831) * t678 + 0.4e1 * t729 * (0.2e1 * t986 + (t812 + t955) * t847 + t841 * t919 - Ifges(3,4)) * t982 + (t897 + (((-qJ(2,3) / 0.4e1 + t869) * mrSges(3,1) + t860) * t729 + ((t813 + t989) * t1010 + t920) * t711) * t847 + t913 / 0.2e1 - t841 * (-t765 * t729 + 0.2e1 * t926) / 0.4e1) * t729 * t1009 - 0.4e1 * t729 * t916 - 0.2e1 * ((((-qJ(2,3) / 0.2e1 + t870) * mrSges(3,2) + t859) * t729 + t926) * t831 + t705) * t729 * t847 - (t711 * t934 - t768 * t729) * t729 * t961 + 0.2e1 * t711 * (Ifges(3,4) * t729 + t756 * t687);
t683 = t710 - t693 / 0.2e1 - t692 / 0.2e1 - t721 / 0.2e1;
t688 = t693 + t692 + t721;
t779 = 0.2e1 * t821;
t874 = qJ(2,2) ^ 2;
t670 = ((t683 * t1015 + ((qJ(2,2) * t1008) - t876 * cos(t779) - (2 * t874) + t880) * t712 / 0.2e1 + (t751 + t1014) * t688) * t712 + ((t730 * t817 * t785 - 0.2e1 * t683 * t788 + (-cos(t868 + qJ(3,2)) * pkin(2) - t997) * t712) * t712 - (-t712 * t817 * sin(t779) / 0.2e1 + t730 * t751) * t782 * t730) * pkin(3)) * t799;
t679 = (-pkin(3) * t975 + (-t710 + 0.2e1 * t688 + (-t790 - t1003) * t712) * t712) * t799;
t706 = t712 * t954;
t715 = -(-t766 * t849 - t769 * t843) * t832 - t831 * (t843 * t766 - t769 * t849);
t892 = mrSges(3,1) * t849 - t991;
t939 = t796 + t793;
t718 = -(-t892 - t983) * t832 - (mrSges(2,2) + t939) * t831 + t814;
t757 = t858 * qJ(2,2) + t917;
t896 = t957 * t979;
t909 = t785 * t975;
t922 = pkin(2) * t991;
t912 = (t922 - t836) * t979;
t915 = t712 * t985;
t925 = mrSges(3,1) * t710;
t929 = 0.4e1 * t988;
t942 = (-((t824 + t929) * t849 - 0.2e1 * t922 + t894) * t819 - (t828 * t931 + t898 * t849 + (t983 - t991) * t1013) * t832 - t957 + (t793 + mrSges(2,2)) * t933 - (t834 ^ 2) * m(3) - m(2) * t874 + t884) * t679 + t718 * t670 - t715 * t909 + 0.2e1 * (Ifges(3,4) * t949 + t819 * t957 - t837 * qJ(2,2) + (-(-pkin(2) * t793 + t836 * t949 + t995) * t832 + pkin(1) * t796) * t831) * t679 + 0.4e1 * t730 * (0.2e1 * t985 + (t812 + t954) * t849 + t843 * t919 - Ifges(3,4)) * t980 + (t896 + (((-qJ(2,2) / 0.4e1 + t869) * mrSges(3,1) + t860) * t730 + ((t813 + t988) * t1010 + t920) * t712) * t849 + t912 / 0.2e1 - t843 * (-t766 * t730 + 0.2e1 * t925) / 0.4e1) * t730 * t1009 - 0.4e1 * t730 * t915 - 0.2e1 * ((((-qJ(2,2) / 0.2e1 + t870) * mrSges(3,2) + t859) * t730 + t925) * t831 + t706) * t730 * t849 - (t712 * t934 - t769 * t730) * t730 * t960 + 0.2e1 * t712 * (Ifges(3,4) * t730 + t757 * t688);
t681 = t708 - t695 / 0.2e1 - t694 / 0.2e1 - t722 / 0.2e1;
t689 = t695 + t694 + t722;
t780 = 0.2e1 * t822;
t875 = qJ(2,1) ^ 2;
t671 = ((t681 * t1015 + ((qJ(2,1) * t1008) - t876 * cos(t780) - (2 * t875) + t880) * t713 / 0.2e1 + (t752 + t1014) * t689) * t713 + ((t731 * t818 * t786 - 0.2e1 * t681 * t789 + (-cos(t868 + qJ(3,1)) * pkin(2) - t996) * t713) * t713 - (-t713 * t818 * sin(t780) / 0.2e1 + t731 * t752) * t783 * t731) * pkin(3)) * t800;
t680 = (-pkin(3) * t974 + (-t708 + 0.2e1 * t689 + (-t790 - t1002) * t713) * t713) * t800;
t707 = t713 * t953;
t716 = -(-t767 * t851 - t770 * t845) * t832 - t831 * (t845 * t767 - t770 * t851);
t891 = mrSges(3,1) * t851 - t990;
t938 = t797 + t794;
t719 = -(-t891 - t983) * t832 - (mrSges(2,2) + t938) * t831 + t814;
t758 = t858 * qJ(2,1) + t917;
t895 = t956 * t977;
t908 = t786 * t974;
t921 = pkin(2) * t990;
t911 = (t921 - t836) * t977;
t914 = t713 * t984;
t924 = mrSges(3,1) * t708;
t928 = 0.4e1 * t987;
t941 = (-((t824 + t928) * t851 - 0.2e1 * t921 + t894) * t819 - (t829 * t931 + t898 * t851 + (t983 - t990) * t1013) * t832 - t956 + (t794 + mrSges(2,2)) * t933 - (t835 ^ 2) * m(3) - m(2) * t875 + t884) * t680 + t719 * t671 - t716 * t908 + 0.2e1 * (Ifges(3,4) * t947 + t819 * t956 - t837 * qJ(2,1) + (-(-pkin(2) * t794 + t836 * t947 + t995) * t832 + pkin(1) * t797) * t831) * t680 + 0.4e1 * t731 * (0.2e1 * t984 + (t812 + t953) * t851 + t845 * t919 - Ifges(3,4)) * t978 + (t895 + (((-qJ(2,1) / 0.4e1 + t869) * mrSges(3,1) + t860) * t731 + ((t813 + t987) * t1010 + t920) * t713) * t851 + t911 / 0.2e1 - t845 * (-t767 * t731 + 0.2e1 * t924) / 0.4e1) * t731 * t1009 - 0.4e1 * t731 * t914 - 0.2e1 * ((((-qJ(2,1) / 0.2e1 + t870) * mrSges(3,2) + t859) * t731 + t924) * t831 + t707) * t731 * t851 - (t713 * t934 - t770 * t731) * t731 * t959 + 0.2e1 * t713 * (Ifges(3,4) * t731 + t758 * t689);
t918 = -t1006 / 0.4e1;
t904 = t943 * t781;
t903 = t942 * t782;
t902 = t941 * t783;
t901 = -t858 * t669 + t717 * t678 - (t756 * t711 + 0.2e1 * (-t893 * t831 - t940 * t832) * t729) * t711;
t900 = -t858 * t670 + t718 * t679 - (t757 * t712 + 0.2e1 * (-t892 * t831 - t939 * t832) * t730) * t712;
t899 = -t858 * t671 + t719 * t680 - (t758 * t713 + 0.2e1 * (-t891 * t831 - t938 * t832) * t731) * t713;
t807 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t811 = -t1005 / 0.4e1;
t861 = -Ifges(3,4) / 0.2e1;
t890 = t781 * ((-0.4e1 * (t986 + (t807 * t841 + t811) * t847 + t841 * t918 + t861) * t982 + (-0.2e1 * t897 + (mrSges(3,2) * t685 + (t930 + t1006) * t981) * t847 - t913 + mrSges(3,1) * t952) * t832 + 0.2e1 * t916 + (t685 * t994 + t705) * t847 - t952 * t993 - Ifges(3,4) * t711) * t711 + Ifges(3,3) * t910 + t714 * t678);
t889 = t782 * ((-0.4e1 * (t985 + (t807 * t843 + t811) * t849 + t843 * t918 + t861) * t980 + (-0.2e1 * t896 + (mrSges(3,2) * t686 + (t929 + t1006) * t979) * t849 - t912 + mrSges(3,1) * t950) * t832 + 0.2e1 * t915 + (t686 * t994 + t706) * t849 - t950 * t993 - Ifges(3,4) * t712) * t712 + Ifges(3,3) * t909 + t715 * t679);
t888 = t783 * ((-0.4e1 * (t984 + (t807 * t845 + t811) * t851 + t845 * t918 + t861) * t978 + (-0.2e1 * t895 + (mrSges(3,2) * t684 + (t928 + t1006) * t977) * t851 - t911 + mrSges(3,1) * t948) * t832 + 0.2e1 * t914 + (t684 * t994 + t707) * t851 - t948 * t993 - Ifges(3,4) * t713) * t713 + Ifges(3,3) * t908 + t716 * t680);
t887 = t901 * t747;
t886 = t900 * t748;
t885 = t899 * t749;
t1 = [(t704 * t885 + t743 * t902) * t800 + (t702 * t886 + t741 * t903) * t799 + (t700 * t887 + t739 * t904) * t798 + (t801 * t890 + t802 * t889 + t803 * t888) * t877; (t703 * t885 + t742 * t902) * t800 + (t701 * t886 + t740 * t903) * t799 + (t699 * t887 + t738 * t904) * t798 + (t804 * t890 + t805 * t889 + t806 * t888) * t877; (t899 * t734 - t941 * t846) * t800 + (t900 * t733 - t942 * t844) * t799 + (t901 * t732 - t943 * t842) * t798;];
taucX  = t1;
