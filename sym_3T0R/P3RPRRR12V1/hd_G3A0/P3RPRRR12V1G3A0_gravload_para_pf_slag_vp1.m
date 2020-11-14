% Calculate Gravitation load for parallel robot
% P3RPRRR12V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:02
% EndTime: 2020-08-06 18:28:02
% DurationCPUTime: 0.41s
% Computational Cost: add. (423->112), mult. (648->199), div. (36->7), fcn. (378->18), ass. (0->71)
t738 = pkin(1) + pkin(5);
t737 = m(3) / pkin(3);
t695 = legFrame(3,2);
t685 = sin(t695);
t688 = cos(t695);
t673 = t688 * g(1) - t685 * g(2);
t699 = sin(qJ(1,3));
t705 = cos(qJ(1,3));
t665 = g(3) * t705 + t673 * t699;
t670 = t685 * g(1) + t688 * g(2);
t698 = sin(qJ(3,3));
t692 = 0.1e1 / t698;
t704 = cos(qJ(3,3));
t736 = ((t665 * rSges(3,1) - rSges(3,2) * t670) * t704 - t698 * (rSges(3,1) * t670 + t665 * rSges(3,2))) * t692;
t696 = legFrame(2,2);
t686 = sin(t696);
t689 = cos(t696);
t674 = t689 * g(1) - t686 * g(2);
t701 = sin(qJ(1,2));
t707 = cos(qJ(1,2));
t666 = g(3) * t707 + t674 * t701;
t671 = t686 * g(1) + t689 * g(2);
t700 = sin(qJ(3,2));
t693 = 0.1e1 / t700;
t706 = cos(qJ(3,2));
t735 = ((t666 * rSges(3,1) - rSges(3,2) * t671) * t706 - t700 * (rSges(3,1) * t671 + t666 * rSges(3,2))) * t693;
t697 = legFrame(1,2);
t687 = sin(t697);
t690 = cos(t697);
t675 = t690 * g(1) - t687 * g(2);
t703 = sin(qJ(1,1));
t709 = cos(qJ(1,1));
t667 = g(3) * t709 + t675 * t703;
t672 = t687 * g(1) + t690 * g(2);
t702 = sin(qJ(3,1));
t694 = 0.1e1 / t702;
t708 = cos(qJ(3,1));
t734 = ((t667 * rSges(3,1) - rSges(3,2) * t672) * t708 - t702 * (rSges(3,1) * t672 + t667 * rSges(3,2))) * t694;
t682 = t698 * pkin(3) + qJ(2,3);
t679 = 0.1e1 / t682;
t733 = t665 * t679;
t683 = t700 * pkin(3) + qJ(2,2);
t680 = 0.1e1 / t683;
t732 = t666 * t680;
t684 = t702 * pkin(3) + qJ(2,1);
t681 = 0.1e1 / t684;
t731 = t667 * t681;
t669 = (rSges(3,3) + t738) * m(3) + (pkin(1) - rSges(2,2)) * m(2) + m(1) * rSges(1,1);
t668 = g(3) * t669;
t710 = m(1) * rSges(1,2);
t676 = -qJ(2,3) * m(3) + (-rSges(2,3) - qJ(2,3)) * m(2) + t710;
t730 = t679 * ((t676 * t673 + t668) * t705 + (-g(3) * t676 + t669 * t673) * t699 + (g(3) * t699 - t673 * t705) * m(3) * (t698 * rSges(3,1) + rSges(3,2) * t704));
t677 = -qJ(2,2) * m(3) + (-rSges(2,3) - qJ(2,2)) * m(2) + t710;
t729 = t680 * ((t677 * t674 + t668) * t707 + (-g(3) * t677 + t669 * t674) * t701 + (g(3) * t701 - t674 * t707) * m(3) * (t700 * rSges(3,1) + rSges(3,2) * t706));
t678 = -qJ(2,1) * m(3) + (-rSges(2,3) - qJ(2,1)) * m(2) + t710;
t728 = t681 * ((t678 * t675 + t668) * t709 + (-g(3) * t678 + t669 * t675) * t703 + (g(3) * t703 - t675 * t709) * m(3) * (t702 * rSges(3,1) + rSges(3,2) * t708));
t727 = t704 * qJ(2,3);
t726 = t706 * qJ(2,2);
t725 = t708 * qJ(2,1);
t724 = t692 * t733;
t723 = t693 * t732;
t722 = t694 * t731;
t721 = t705 * t730;
t720 = t707 * t729;
t719 = t709 * t728;
t691 = pkin(6) + t738;
t715 = qJ(2,1) * t703 + t691 * t709;
t714 = qJ(2,2) * t701 + t691 * t707;
t713 = qJ(2,3) * t699 + t691 * t705;
t711 = -m(2) - m(3);
t1 = [-m(4) * g(1) + t688 * t721 + t689 * t720 + t690 * t719 + (t685 * t736 + t686 * t735 + t687 * t734) * t737 + ((t715 * t690 * t702 + t687 * t725 + (t687 * t708 * t702 + (-t708 ^ 2 + 0.1e1) * t690 * t703) * pkin(3)) * t722 + (t714 * t689 * t700 + t686 * t726 + (t686 * t706 * t700 + (-t706 ^ 2 + 0.1e1) * t689 * t701) * pkin(3)) * t723 + (t713 * t688 * t698 + t685 * t727 + (t685 * t704 * t698 + (-t704 ^ 2 + 0.1e1) * t688 * t699) * pkin(3)) * t724) * t711; -t685 * t721 - t686 * t720 - t687 * t719 - m(4) * g(2) + (((t690 * pkin(3) * t708 - t715 * t687) * t702 + t703 * pkin(3) * (t708 - 0.1e1) * (t708 + 0.1e1) * t687 + t690 * t725) * t722 + ((t689 * pkin(3) * t706 - t714 * t686) * t700 + t701 * pkin(3) * (t706 - 0.1e1) * (t706 + 0.1e1) * t686 + t689 * t726) * t723 + ((t688 * pkin(3) * t704 - t713 * t685) * t698 + t699 * pkin(3) * (t704 - 0.1e1) * (t704 + 0.1e1) * t685 + t688 * t727) * t724) * t711 + (t688 * t736 + t689 * t735 + t690 * t734) * t737; -t699 * t730 - t701 * t729 - t703 * t728 - m(4) * g(3) + ((t684 * t709 - t691 * t703) * t731 + (t683 * t707 - t691 * t701) * t732 + (t682 * t705 - t691 * t699) * t733) * t711;];
taugX  = t1;
