% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR1G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:13
% EndTime: 2020-03-09 20:34:14
% DurationCPUTime: 0.36s
% Computational Cost: add. (161->65), mult. (367->146), div. (108->10), fcn. (450->18), ass. (0->77)
t833 = legFrame(1,3);
t815 = sin(t833);
t818 = cos(t833);
t812 = g(1) * t818 + g(2) * t815;
t839 = sin(qJ(2,1));
t845 = cos(qJ(2,1));
t799 = -g(3) * t845 + t812 * t839;
t832 = legFrame(2,3);
t814 = sin(t832);
t817 = cos(t832);
t811 = g(1) * t817 + g(2) * t814;
t837 = sin(qJ(2,2));
t843 = cos(qJ(2,2));
t797 = -g(3) * t843 + t811 * t837;
t831 = legFrame(3,3);
t813 = sin(t831);
t816 = cos(t831);
t810 = g(1) * t816 + g(2) * t813;
t835 = sin(qJ(2,3));
t841 = cos(qJ(2,3));
t795 = -g(3) * t841 + t810 * t835;
t807 = g(1) * t813 - g(2) * t816;
t840 = cos(qJ(3,3));
t825 = 0.1e1 / t840;
t834 = sin(qJ(3,3));
t852 = g(3) * t835 + t810 * t841;
t879 = (t807 * t834 + t840 * t852) * t825;
t808 = g(1) * t814 - g(2) * t817;
t842 = cos(qJ(3,2));
t827 = 0.1e1 / t842;
t836 = sin(qJ(3,2));
t851 = g(3) * t837 + t811 * t843;
t878 = (t808 * t836 + t842 * t851) * t827;
t809 = g(1) * t815 - g(2) * t818;
t844 = cos(qJ(3,1));
t829 = 0.1e1 / t844;
t838 = sin(qJ(3,1));
t850 = g(3) * t839 + t812 * t845;
t877 = (t809 * t838 + t844 * t850) * t829;
t822 = 0.1e1 / t835;
t876 = t795 * t822;
t823 = 0.1e1 / t837;
t875 = t797 * t823;
t824 = 0.1e1 / t839;
t874 = t799 * t824;
t870 = t822 * t825;
t869 = t822 / t840 ^ 2;
t868 = t823 * t827;
t867 = t823 / t842 ^ 2;
t866 = t824 * t829;
t865 = t824 / t844 ^ 2;
t864 = t834 * t841;
t863 = t836 * t843;
t862 = t838 * t845;
t861 = t840 * t841;
t860 = t842 * t843;
t859 = t844 * t845;
t801 = -t813 * t864 - t816 * t840;
t858 = t801 * t869;
t802 = -t814 * t863 - t817 * t842;
t857 = t802 * t867;
t803 = -t815 * t862 - t818 * t844;
t856 = t803 * t865;
t804 = -t813 * t840 + t816 * t864;
t855 = t804 * t869;
t805 = -t814 * t842 + t817 * t863;
t854 = t805 * t867;
t806 = -t815 * t844 + t818 * t862;
t853 = t806 * t865;
t849 = t795 * t834 * t869;
t848 = t797 * t836 * t867;
t847 = t799 * t838 * t865;
t846 = 0.1e1 / pkin(2);
t787 = -t809 * t844 + t838 * t850;
t785 = -t808 * t842 + t836 * t851;
t783 = -t807 * t840 + t834 * t852;
t1 = [(-(t815 * t838 + t818 * t859) * t866 - (t814 * t836 + t817 * t860) * t868 - (t813 * t834 + t816 * t861) * t870) * g(3), 0, (t795 * t858 + t797 * t857 + t799 * t856) * t846, (t850 * t856 + t851 * t857 + t852 * t858) * t846, 0, 0, 0, 0, 0, ((t787 * t815 + t803 * t874) * t829 + (t785 * t814 + t802 * t875) * t827 + (t783 * t813 + t801 * t876) * t825) * t846, (-t801 * t849 - t802 * t848 - t803 * t847 + t813 * t879 + t814 * t878 + t815 * t877) * t846, -g(1); (-(t815 * t859 - t818 * t838) * t866 - (t814 * t860 - t817 * t836) * t868 - (t813 * t861 - t816 * t834) * t870) * g(3), 0, (t795 * t855 + t797 * t854 + t799 * t853) * t846, (t850 * t853 + t851 * t854 + t852 * t855) * t846, 0, 0, 0, 0, 0, ((-t787 * t818 + t806 * t874) * t829 + (-t785 * t817 + t805 * t875) * t827 + (-t783 * t816 + t804 * t876) * t825) * t846, (-t804 * t849 - t805 * t848 - t806 * t847 - t816 * t879 - t817 * t878 - t818 * t877) * t846, -g(2); -0.3e1 * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
tau_reg  = t1;
