% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G4A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:25:14
% EndTime: 2020-08-06 17:25:21
% DurationCPUTime: 6.56s
% Computational Cost: add. (36948->282), mult. (97179->604), div. (6885->8), fcn. (124392->34), ass. (0->264)
t867 = cos(qJ(3,1));
t841 = 0.1e1 / t867;
t868 = cos(qJ(2,1));
t862 = sin(qJ(2,1));
t908 = t862 * t867;
t799 = pkin(2) * t908 - t868 * pkin(5);
t843 = sin(pkin(3));
t861 = sin(qJ(3,1));
t845 = cos(pkin(3));
t965 = pkin(2) * t845;
t778 = t799 * t843 + t861 * t965;
t978 = 0.1e1 / t778;
t946 = t978 * t841;
t865 = cos(qJ(3,2));
t839 = 0.1e1 / t865;
t866 = cos(qJ(2,2));
t860 = sin(qJ(2,2));
t910 = t860 * t865;
t798 = pkin(2) * t910 - t866 * pkin(5);
t859 = sin(qJ(3,2));
t777 = t798 * t843 + t859 * t965;
t979 = 0.1e1 / t777;
t947 = t979 * t839;
t857 = sin(qJ(3,3));
t863 = cos(qJ(3,3));
t895 = t863 * mrSges(3,1) - t857 * mrSges(3,2);
t894 = t865 * mrSges(3,1) - t859 * mrSges(3,2);
t893 = t867 * mrSges(3,1) - t861 * mrSges(3,2);
t981 = -2 * Ifges(3,4);
t858 = sin(qJ(2,3));
t912 = t858 * t863;
t898 = t843 * t912;
t921 = t845 * t857;
t864 = cos(qJ(2,3));
t926 = t843 * t864;
t980 = 0.1e1 / (-pkin(5) * t926 + (t898 + t921) * pkin(2));
t836 = t863 ^ 2;
t977 = 0.2e1 * t836;
t838 = t865 ^ 2;
t976 = 0.2e1 * t838;
t840 = t867 ^ 2;
t975 = 0.2e1 * t840;
t974 = Ifges(3,5) / 0.2e1;
t973 = -Ifges(3,6) / 0.2e1;
t853 = mrSges(2,2) - mrSges(3,3);
t972 = t853 / 0.2e1;
t846 = legFrame(3,3);
t814 = sin(t846);
t820 = cos(t846);
t842 = sin(pkin(6));
t844 = cos(pkin(6));
t782 = t814 * t844 + t820 * t842;
t785 = -t814 * t842 + t820 * t844;
t849 = legFrame(3,1);
t817 = sin(t849);
t823 = cos(t849);
t854 = legFrame(3,2);
t829 = sin(t854);
t939 = t823 * t829;
t744 = t782 * t939 + t817 * t785;
t745 = -t817 * t782 + t785 * t939;
t917 = t845 * t864;
t920 = t845 * t858;
t964 = pkin(2) * t863;
t717 = (-t858 * t744 + t745 * t917) * t964 + pkin(5) * (t864 * t744 + t745 * t920);
t942 = t817 * t829;
t750 = t782 * t823 + t785 * t942;
t753 = t782 * t942 - t785 * t823;
t718 = -(t750 * t917 - t753 * t858) * t964 - pkin(5) * (t750 * t920 + t864 * t753);
t738 = (-t858 * t782 + t785 * t917) * t964 + pkin(5) * (t864 * t782 + t785 * t920);
t869 = xDP(3);
t870 = xDP(2);
t874 = 0.1e1 / pkin(2);
t832 = cos(t854);
t871 = xDP(1);
t935 = t832 * t871;
t837 = 0.1e1 / t863;
t948 = t980 * t837;
t702 = (t717 * t869 + t718 * t870 - t738 * t935) * t874 * t948;
t971 = pkin(2) * t702;
t847 = legFrame(2,3);
t815 = sin(t847);
t821 = cos(t847);
t783 = t815 * t844 + t821 * t842;
t786 = -t815 * t842 + t821 * t844;
t850 = legFrame(2,1);
t818 = sin(t850);
t824 = cos(t850);
t855 = legFrame(2,2);
t830 = sin(t855);
t938 = t824 * t830;
t746 = t783 * t938 + t818 * t786;
t747 = -t818 * t783 + t786 * t938;
t916 = t845 * t866;
t919 = t845 * t860;
t963 = pkin(2) * t865;
t719 = (-t860 * t746 + t747 * t916) * t963 + pkin(5) * (t866 * t746 + t747 * t919);
t941 = t818 * t830;
t751 = t783 * t824 + t786 * t941;
t754 = t783 * t941 - t786 * t824;
t720 = -(t751 * t916 - t754 * t860) * t963 - pkin(5) * (t751 * t919 + t866 * t754);
t739 = (-t860 * t783 + t786 * t916) * t963 + pkin(5) * (t866 * t783 + t786 * t919);
t833 = cos(t855);
t933 = t833 * t871;
t703 = (t719 * t869 + t720 * t870 - t739 * t933) * t874 * t947;
t970 = pkin(2) * t703;
t848 = legFrame(1,3);
t816 = sin(t848);
t822 = cos(t848);
t784 = t816 * t844 + t822 * t842;
t787 = -t816 * t842 + t822 * t844;
t851 = legFrame(1,1);
t819 = sin(t851);
t825 = cos(t851);
t856 = legFrame(1,2);
t831 = sin(t856);
t937 = t825 * t831;
t748 = t784 * t937 + t819 * t787;
t749 = -t819 * t784 + t787 * t937;
t915 = t845 * t868;
t918 = t845 * t862;
t962 = pkin(2) * t867;
t721 = (-t862 * t748 + t749 * t915) * t962 + pkin(5) * (t868 * t748 + t749 * t918);
t940 = t819 * t831;
t752 = t784 * t825 + t787 * t940;
t755 = t784 * t940 - t787 * t825;
t722 = -(t752 * t915 - t755 * t862) * t962 - pkin(5) * (t752 * t918 + t868 * t755);
t740 = (-t862 * t784 + t787 * t915) * t962 + pkin(5) * (t868 * t784 + t787 * t918);
t834 = cos(t856);
t931 = t834 * t871;
t704 = (t721 * t869 + t722 * t870 - t740 * t931) * t874 * t946;
t969 = pkin(2) * t704;
t968 = pkin(2) * t836;
t967 = pkin(2) * t838;
t966 = pkin(2) * t840;
t961 = -Ifges(3,1) - Ifges(2,3);
t789 = t842 * t866 + t844 * t919;
t792 = -t842 * t919 + t844 * t866;
t763 = t789 * t821 + t815 * t792;
t882 = t815 * t789 - t792 * t821;
t925 = t843 * t865;
t724 = (t763 * t938 - t818 * t882) * t859 + t747 * t925;
t727 = (-t763 * t941 - t824 * t882) * t859 - t751 * t925;
t742 = t763 * t859 + t786 * t925;
t709 = (t724 * t869 + t727 * t870 - t742 * t933) * t947;
t872 = pkin(5) ^ 2;
t873 = pkin(2) ^ 2;
t954 = t703 * t859;
t903 = pkin(2) * t954;
t957 = (-pkin(5) * t903 + (t838 * t873 + t872) * t709) * t709;
t790 = t842 * t868 + t844 * t918;
t793 = -t842 * t918 + t844 * t868;
t764 = t790 * t822 + t816 * t793;
t881 = t816 * t790 - t793 * t822;
t923 = t843 * t867;
t725 = (t764 * t937 - t819 * t881) * t861 + t749 * t923;
t728 = (-t764 * t940 - t825 * t881) * t861 - t752 * t923;
t743 = t764 * t861 + t787 * t923;
t710 = (t725 * t869 + t728 * t870 - t743 * t931) * t946;
t953 = t704 * t861;
t902 = pkin(2) * t953;
t956 = (-pkin(5) * t902 + (t840 * t873 + t872) * t710) * t710;
t955 = t702 * t857;
t788 = t842 * t864 + t844 * t920;
t791 = -t842 * t920 + t844 * t864;
t762 = t788 * t820 + t814 * t791;
t883 = t814 * t788 - t791 * t820;
t927 = t843 * t863;
t723 = (t762 * t939 - t817 * t883) * t857 + t745 * t927;
t726 = (-t762 * t942 - t823 * t883) * t857 - t750 * t927;
t741 = t762 * t857 + t785 * t927;
t797 = pkin(2) * t912 - t864 * pkin(5);
t776 = pkin(2) * t921 + t797 * t843;
t773 = 0.1e1 / t776;
t708 = (t723 * t869 + t726 * t870 - t741 * t935) * t837 * t773;
t952 = t708 * t980;
t951 = t708 * t857;
t950 = t709 * t859;
t949 = t710 * t861;
t945 = ((mrSges(2,1) + t895) * t864 - t858 * t853) * t843;
t944 = ((mrSges(2,1) + t894) * t866 - t860 * t853) * t843;
t943 = ((mrSges(2,1) + t893) * t868 - t862 * t853) * t843;
t936 = t832 * t837;
t934 = t833 * t839;
t932 = t834 * t841;
t930 = t843 * t857;
t929 = t843 * t859;
t928 = t843 * t861;
t924 = t843 * t866;
t922 = t843 * t868;
t914 = t845 * t874;
t913 = t857 * t863;
t911 = t859 * t865;
t909 = t861 * t867;
t901 = pkin(5) * t951;
t693 = t901 - t971;
t904 = pkin(2) * t955;
t675 = (((t845 * t702 + t708 * t926) * t968 - (-pkin(5) * t708 + t904) * t898 + t845 * t693) * t952 + (t702 * t926 + (t836 * t845 - t857 * t898 - t845) * t708) * t980 * t971) * t837;
t690 = -pkin(5) * t904 + (t836 * t873 + t872) * t708;
t681 = (t690 * t914 * t952 + (-t702 * t797 * t930 + t845 * (t702 * t968 - t901)) * t773 * t702) * t837;
t684 = (t690 * t708 - t693 * t971) * t980;
t699 = t702 ^ 2;
t705 = t708 ^ 2;
t803 = t857 * mrSges(3,1) + t863 * mrSges(3,2);
t765 = t803 * t843 * t858 - t895 * t845;
t835 = -m(1) - m(2) - m(3);
t907 = -t675 * t945 + t765 * t681 - t835 * t684 + ((-t705 * mrSges(2,1) - t895 * (t705 + t699)) * t858 - 0.2e1 * t708 * t864 * (t803 * t702 + t708 * t972)) * t843 - t845 * t699 * t803;
t900 = pkin(5) * t950;
t694 = t900 - t970;
t897 = t843 * t910;
t676 = (((t845 * t703 + t709 * t924) * t967 - (-pkin(5) * t709 + t903) * t897 + t845 * t694) * t709 + (t703 * t924 + (t838 * t845 - t859 * t897 - t845) * t709) * t970) * t947;
t682 = (t914 * t957 + (-t703 * t798 * t929 + t845 * (t703 * t967 - t900)) * t703) * t947;
t685 = (t694 * t970 - t957) * t979;
t700 = t703 ^ 2;
t706 = t709 ^ 2;
t804 = t859 * mrSges(3,1) + t865 * mrSges(3,2);
t766 = t804 * t843 * t860 - t894 * t845;
t906 = -t676 * t944 + t766 * t682 + t835 * t685 + ((-t706 * mrSges(2,1) - t894 * (t706 + t700)) * t860 - 0.2e1 * t709 * t866 * (t804 * t703 + t709 * t972)) * t843 - t845 * t700 * t804;
t899 = pkin(5) * t949;
t695 = t899 - t969;
t896 = t843 * t908;
t677 = (((t845 * t704 + t710 * t922) * t966 - (-pkin(5) * t710 + t902) * t896 + t845 * t695) * t710 - (-t704 * t922 + (-t840 * t845 + t861 * t896 + t845) * t710) * t969) * t946;
t683 = (t914 * t956 + (-t704 * t799 * t928 + t845 * (t704 * t966 - t899)) * t704) * t946;
t686 = (t695 * t969 - t956) * t978;
t701 = t704 ^ 2;
t707 = t710 ^ 2;
t805 = t861 * mrSges(3,1) + t867 * mrSges(3,2);
t767 = t805 * t843 * t862 - t893 * t845;
t905 = -t677 * t943 + t767 * t683 + t835 * t686 + ((-t707 * mrSges(2,1) - t893 * (t707 + t701)) * t862 - 0.2e1 * t710 * t868 * (t805 * t704 + t710 * t972)) * t843 - t845 * t701 * t805;
t806 = -Ifges(3,5) * t857 - Ifges(3,6) * t863;
t852 = Ifges(3,1) - Ifges(3,2);
t892 = 0.2e1 * ((t702 * t974 + t852 * t951) * t863 + t955 * t973 + (t977 - 0.1e1) * t708 * Ifges(3,4)) * t702 + t684 * t945 + (t852 * t836 + t913 * t981 + t961) * t675 + t806 * t681;
t807 = -Ifges(3,5) * t859 - Ifges(3,6) * t865;
t891 = 0.2e1 * ((t703 * t974 + t852 * t950) * t865 + t954 * t973 + (t976 - 0.1e1) * t709 * Ifges(3,4)) * t703 - t685 * t944 + (t852 * t838 + t911 * t981 + t961) * t676 + t807 * t682;
t808 = -Ifges(3,5) * t861 - Ifges(3,6) * t867;
t890 = 0.2e1 * ((t704 * t974 + t852 * t949) * t867 + t953 * t973 + (t975 - 0.1e1) * t710 * Ifges(3,4)) * t704 - t686 * t943 + (t852 * t840 + t909 * t981 + t961) * t677 + t808 * t683;
t889 = Ifges(3,3) * t681 - t806 * t675 + t765 * t684 + t705 * (Ifges(3,4) * t977 + t852 * t913 - Ifges(3,4));
t888 = Ifges(3,3) * t682 - t807 * t676 - t766 * t685 + t706 * (Ifges(3,4) * t976 + t852 * t911 - Ifges(3,4));
t887 = Ifges(3,3) * t683 - t808 * t677 - t767 * t686 + t707 * (Ifges(3,4) * t975 + t852 * t909 - Ifges(3,4));
t886 = t892 * t837;
t885 = t891 * t839;
t884 = t890 * t841;
t880 = pkin(2) * t930 - t797 * t845;
t879 = pkin(2) * t929 - t798 * t845;
t878 = pkin(2) * t928 - t799 * t845;
t877 = t889 * t948;
t876 = t888 * t947;
t875 = t887 * t946;
t802 = pkin(5) * t862 + t868 * t962;
t801 = pkin(5) * t860 + t866 * t963;
t800 = pkin(5) * t858 + t864 * t964;
t761 = -t842 * t802 + t844 * t878;
t760 = -t842 * t801 + t844 * t879;
t759 = -t842 * t800 + t844 * t880;
t758 = t802 * t844 + t842 * t878;
t757 = t801 * t844 + t842 * t879;
t756 = t800 * t844 + t842 * t880;
t737 = -t758 * t816 + t761 * t822;
t736 = -t757 * t815 + t760 * t821;
t735 = -t756 * t814 + t759 * t820;
t731 = -t834 * t778 + (t758 * t822 + t761 * t816) * t831;
t730 = -t833 * t777 + (t757 * t821 + t760 * t815) * t830;
t729 = -t832 * t776 + (t756 * t820 + t759 * t814) * t829;
t1 = [(-t890 * t743 * t932 + t905 * ((t784 * t878 + t787 * t802) * t834 + t831 * t778)) * t978 + (-t891 * t742 * t934 + t906 * ((t783 * t879 + t786 * t801) * t833 + t830 * t777)) * t979 + (-t892 * t741 * t936 + t907 * ((t782 * t880 + t785 * t800) * t832 + t829 * t776)) * t773 + (t738 * t889 * t936 * t980 + t739 * t888 * t934 * t979 + t740 * t887 * t932 * t978) * t874; (t728 * t884 + t905 * (t731 * t819 - t825 * t737)) * t978 + (t727 * t885 + t906 * (t730 * t818 - t824 * t736)) * t979 + (t726 * t886 + t907 * (t729 * t817 - t823 * t735)) * t773 + (-t718 * t877 - t720 * t876 - t722 * t875) * t874; (t725 * t884 + t905 * (-t731 * t825 - t819 * t737)) * t978 + (t724 * t885 + t906 * (-t730 * t824 - t818 * t736)) * t979 + (t723 * t886 + t907 * (-t729 * t823 - t817 * t735)) * t773 + (-t717 * t877 - t719 * t876 - t721 * t875) * t874;];
taucX  = t1;
