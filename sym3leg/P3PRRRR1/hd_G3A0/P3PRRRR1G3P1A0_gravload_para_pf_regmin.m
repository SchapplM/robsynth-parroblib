% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR1G3P1A0
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
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G3P1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:28
% EndTime: 2020-03-09 21:02:29
% DurationCPUTime: 0.37s
% Computational Cost: add. (186->72), mult. (369->165), div. (126->10), fcn. (447->18), ass. (0->74)
t735 = legFrame(3,2);
t720 = sin(t735);
t723 = cos(t735);
t714 = t720 * g(1) + t723 * g(2);
t717 = t723 * g(1) - t720 * g(2);
t739 = sin(qJ(2,3));
t745 = cos(qJ(2,3));
t704 = t714 * t745 - t717 * t739;
t726 = 0.1e1 / t739;
t780 = t704 * t726;
t736 = legFrame(2,2);
t721 = sin(t736);
t724 = cos(t736);
t715 = t721 * g(1) + t724 * g(2);
t718 = t724 * g(1) - t721 * g(2);
t741 = sin(qJ(2,2));
t747 = cos(qJ(2,2));
t707 = t715 * t747 - t718 * t741;
t727 = 0.1e1 / t741;
t779 = t707 * t727;
t737 = legFrame(1,2);
t722 = sin(t737);
t725 = cos(t737);
t716 = t722 * g(1) + t725 * g(2);
t719 = t725 * g(1) - t722 * g(2);
t743 = sin(qJ(2,1));
t749 = cos(qJ(2,1));
t710 = t716 * t749 - t719 * t743;
t728 = 0.1e1 / t743;
t778 = t710 * t728;
t777 = t714 * t726;
t776 = t715 * t727;
t775 = t716 * t728;
t744 = cos(qJ(3,3));
t729 = 0.1e1 / t744;
t774 = t726 * t729;
t773 = t726 * t745;
t746 = cos(qJ(3,2));
t731 = 0.1e1 / t746;
t772 = t727 * t731;
t771 = t727 * t747;
t748 = cos(qJ(3,1));
t733 = 0.1e1 / t748;
t770 = t728 * t733;
t769 = t728 * t749;
t768 = t720 * t774;
t767 = t721 * t772;
t766 = t722 * t770;
t765 = t723 * t774;
t764 = t724 * t772;
t763 = t725 * t770;
t738 = sin(qJ(3,3));
t762 = t738 * t774;
t761 = 0.1e1 / t744 ^ 2 * t773;
t740 = sin(qJ(3,2));
t760 = t740 * t772;
t759 = 0.1e1 / t746 ^ 2 * t771;
t742 = sin(qJ(3,1));
t758 = t742 * t770;
t757 = 0.1e1 / t748 ^ 2 * t769;
t756 = t704 * t762;
t755 = t707 * t760;
t754 = t710 * t758;
t753 = t738 * t761;
t752 = t740 * t759;
t751 = t742 * t757;
t750 = 0.1e1 / pkin(2);
t713 = (g(1) * t749 + t743 * g(2)) * t725 + t722 * (g(1) * t743 - g(2) * t749);
t712 = (g(1) * t747 + g(2) * t741) * t724 + t721 * (g(1) * t741 - g(2) * t747);
t711 = (g(1) * t745 + g(2) * t739) * t723 + t720 * (g(1) * t739 - g(2) * t745);
t709 = t716 * t743 + t719 * t749;
t706 = t715 * t741 + t718 * t747;
t703 = t714 * t739 + t717 * t745;
t1 = [-(t743 * t722 + t725 * t749) * t775 - (t741 * t721 + t724 * t747) * t776 - (t739 * t720 + t723 * t745) * t777, 0, (t704 * t765 + t707 * t764 + t710 * t763) * t750, (-t703 * t765 - t706 * t764 - t709 * t763) * t750, 0, 0, 0, 0, 0, (t723 * t780 + t724 * t779 + t725 * t778) * t750, (-t723 * t756 - t724 * t755 - t725 * t754) * t750, -g(1); -(-t722 * t749 + t725 * t743) * t775 - (-t721 * t747 + t724 * t741) * t776 - (-t720 * t745 + t723 * t739) * t777, 0, (-t704 * t768 - t707 * t767 - t710 * t766) * t750, (t703 * t768 + t706 * t767 + t709 * t766) * t750, 0, 0, 0, 0, 0, (-t720 * t780 - t721 * t779 - t722 * t778) * t750, (t720 * t756 + t721 * t755 + t722 * t754) * t750, -g(2); -t714 * t762 - t715 * t760 - t716 * t758, 0, (t704 * t753 + t707 * t752 + t710 * t751) * t750, (-t703 * t753 - t706 * t752 - t709 * t751) * t750, 0, 0, 0, 0, 0, ((-g(3) * t748 + (t710 * t769 + t713) * t742) * t733 + (-g(3) * t746 + (t707 * t771 + t712) * t740) * t731 + (-g(3) * t744 + (t704 * t773 + t711) * t738) * t729) * t750, (-t742 ^ 2 * t710 * t757 + t733 * (g(3) * t742 + t713 * t748) - t740 ^ 2 * t707 * t759 + t731 * (g(3) * t740 + t712 * t746) - t738 ^ 2 * t704 * t761 + t729 * (g(3) * t738 + t711 * t744)) * t750, -g(3);];
tau_reg  = t1;
