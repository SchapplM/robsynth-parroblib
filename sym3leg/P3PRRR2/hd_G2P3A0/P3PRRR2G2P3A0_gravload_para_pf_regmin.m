% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR2G2P3A0
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
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR2G2P3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:42
% EndTime: 2020-03-09 21:21:42
% DurationCPUTime: 0.33s
% Computational Cost: add. (307->82), mult. (432->159), div. (126->5), fcn. (492->24), ass. (0->87)
t714 = qJ(2,3) + qJ(3,3);
t699 = sin(t714);
t720 = sin(qJ(3,3));
t711 = 0.1e1 / t720;
t721 = sin(qJ(2,3));
t760 = (t721 * pkin(1) + pkin(2) * t699) * t711;
t715 = qJ(2,2) + qJ(3,2);
t700 = sin(t715);
t722 = sin(qJ(3,2));
t712 = 0.1e1 / t722;
t723 = sin(qJ(2,2));
t759 = (t723 * pkin(1) + pkin(2) * t700) * t712;
t716 = qJ(2,1) + qJ(3,1);
t701 = sin(t716);
t724 = sin(qJ(3,1));
t713 = 0.1e1 / t724;
t725 = sin(qJ(2,1));
t758 = (t725 * pkin(1) + pkin(2) * t701) * t713;
t757 = t699 * t711;
t756 = t700 * t712;
t755 = t701 * t713;
t717 = legFrame(3,2);
t705 = sin(t717);
t754 = t705 * t711;
t718 = legFrame(2,2);
t706 = sin(t718);
t753 = t706 * t712;
t719 = legFrame(1,2);
t707 = sin(t719);
t752 = t707 * t713;
t708 = cos(t717);
t751 = t708 * t711;
t709 = cos(t718);
t750 = t709 * t712;
t710 = cos(t719);
t749 = t710 * t713;
t748 = t721 * t720;
t747 = t723 * t722;
t746 = t725 * t724;
t726 = cos(qJ(3,3));
t727 = cos(qJ(2,3));
t684 = (pkin(2) * t726 + pkin(1)) * t727 - pkin(2) * t748;
t745 = t684 * t754;
t744 = t684 * t751;
t728 = cos(qJ(3,2));
t729 = cos(qJ(2,2));
t685 = (pkin(2) * t728 + pkin(1)) * t729 - pkin(2) * t747;
t743 = t685 * t753;
t742 = t685 * t750;
t730 = cos(qJ(3,1));
t731 = cos(qJ(2,1));
t686 = (pkin(2) * t730 + pkin(1)) * t731 - pkin(2) * t746;
t741 = t686 * t752;
t740 = t686 * t749;
t687 = t727 * t726 - t748;
t739 = t687 * t754;
t738 = t687 * t751;
t688 = t729 * t728 - t747;
t737 = t688 * t753;
t736 = t688 * t750;
t689 = t731 * t730 - t746;
t735 = t689 * t752;
t734 = t689 * t749;
t733 = 0.1e1 / pkin(1);
t732 = 0.1e1 / pkin(2);
t704 = cos(t716);
t703 = cos(t715);
t702 = cos(t714);
t695 = t710 * g(1) - t707 * g(2);
t694 = t709 * g(1) - t706 * g(2);
t693 = t708 * g(1) - t705 * g(2);
t692 = t707 * g(1) + t710 * g(2);
t691 = t706 * g(1) + t709 * g(2);
t690 = t705 * g(1) + t708 * g(2);
t683 = g(3) * t731 + t692 * t725;
t682 = g(3) * t729 + t691 * t723;
t681 = g(3) * t727 + t690 * t721;
t680 = -g(3) * t725 + t692 * t731;
t679 = -g(3) * t723 + t691 * t729;
t678 = -g(3) * t721 + t690 * t727;
t677 = g(3) * t704 + t692 * t701;
t676 = g(3) * t703 + t691 * t700;
t675 = g(3) * t702 + t690 * t699;
t674 = -g(3) * t701 + t692 * t704;
t673 = -g(3) * t700 + t691 * t703;
t672 = -g(3) * t699 + t690 * t702;
t1 = [-t708 * t693 - t709 * t694 - t710 * t695, 0, (t681 * t739 + t682 * t737 + t683 * t735) * t733, (t678 * t739 + t679 * t737 + t680 * t735) * t733, 0, (t675 * t739 + t676 * t737 + t677 * t735 + (-t675 * t745 - t676 * t743 - t677 * t741) * t732) * t733, (t672 * t739 + t673 * t737 + t674 * t735 + (-t672 * t745 - t673 * t743 - t674 * t741) * t732) * t733, -g(1); t705 * t693 + t706 * t694 + t707 * t695, 0, (t681 * t738 + t682 * t736 + t683 * t734) * t733, (t678 * t738 + t679 * t736 + t680 * t734) * t733, 0, (t675 * t738 + t676 * t736 + t677 * t734 + (-t675 * t744 - t676 * t742 - t677 * t740) * t732) * t733, (t672 * t738 + t673 * t736 + t674 * t734 + (-t672 * t744 - t673 * t742 - t674 * t740) * t732) * t733, -g(2); 0, 0, (-t681 * t757 - t682 * t756 - t683 * t755) * t733, (-t678 * t757 - t679 * t756 - t680 * t755) * t733, 0, (-t675 * t757 - t676 * t756 - t677 * t755 + (t675 * t760 + t676 * t759 + t677 * t758) * t732) * t733, (-t672 * t757 - t673 * t756 - t674 * t755 + (t672 * t760 + t673 * t759 + t674 * t758) * t732) * t733, -g(3);];
tau_reg  = t1;
