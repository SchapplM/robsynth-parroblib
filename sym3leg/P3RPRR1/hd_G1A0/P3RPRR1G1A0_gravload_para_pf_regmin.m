% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRR1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:20
% EndTime: 2020-03-09 21:23:20
% DurationCPUTime: 0.28s
% Computational Cost: add. (527->79), mult. (421->115), div. (54->4), fcn. (402->42), ass. (0->79)
t621 = legFrame(3,3);
t650 = cos(t621);
t622 = legFrame(2,3);
t649 = cos(t622);
t623 = legFrame(1,3);
t648 = cos(t623);
t633 = pkin(7) + qJ(3,3);
t588 = 0.1e1 / (pkin(1) * sin(t633) + sin(qJ(3,3)) * pkin(2));
t612 = t621 + qJ(1,3);
t600 = pkin(7) + t612;
t597 = qJ(3,3) + t600;
t591 = sin(t597);
t647 = (-pkin(1) * sin(t612) - pkin(2) * sin(t600) - pkin(3) * t591) * t588;
t634 = pkin(7) + qJ(3,2);
t589 = 0.1e1 / (pkin(1) * sin(t634) + sin(qJ(3,2)) * pkin(2));
t613 = t622 + qJ(1,2);
t601 = pkin(7) + t613;
t598 = qJ(3,2) + t601;
t592 = sin(t598);
t646 = (-pkin(1) * sin(t613) - pkin(2) * sin(t601) - pkin(3) * t592) * t589;
t635 = pkin(7) + qJ(3,1);
t590 = 0.1e1 / (pkin(1) * sin(t635) + sin(qJ(3,1)) * pkin(2));
t614 = t623 + qJ(1,1);
t602 = pkin(7) + t614;
t599 = qJ(3,1) + t602;
t593 = sin(t599);
t645 = (-pkin(1) * sin(t614) - pkin(2) * sin(t602) - pkin(3) * t593) * t590;
t594 = cos(t597);
t644 = (-pkin(1) * cos(t612) - pkin(2) * cos(t600) - pkin(3) * t594) * t588;
t595 = cos(t598);
t643 = (-pkin(1) * cos(t613) - pkin(2) * cos(t601) - pkin(3) * t595) * t589;
t596 = cos(t599);
t642 = (-pkin(1) * cos(t614) - pkin(2) * cos(t602) - pkin(3) * t596) * t590;
t641 = t588 * t591;
t640 = t588 * t594;
t639 = t589 * t592;
t638 = t589 * t595;
t637 = t590 * t593;
t636 = t590 * t596;
t618 = sin(t621);
t582 = t618 * g(1) - t650 * g(2);
t585 = t650 * g(1) + t618 * g(2);
t624 = sin(qJ(1,3));
t627 = cos(qJ(1,3));
t570 = t582 * t627 + t585 * t624;
t619 = sin(t622);
t583 = t619 * g(1) - t649 * g(2);
t586 = t649 * g(1) + t619 * g(2);
t625 = sin(qJ(1,2));
t628 = cos(qJ(1,2));
t571 = t583 * t628 + t586 * t625;
t620 = sin(t623);
t584 = t620 * g(1) - t648 * g(2);
t587 = t648 * g(1) + t620 * g(2);
t626 = sin(qJ(1,1));
t629 = cos(qJ(1,1));
t572 = t584 * t629 + t587 * t626;
t632 = t570 * t641 + t571 * t639 + t572 * t637;
t631 = t570 * t640 + t571 * t638 + t572 * t636;
t630 = 0.1e1 / pkin(3);
t617 = qJ(1,1) + t635;
t616 = qJ(1,2) + t634;
t615 = qJ(1,3) + t633;
t608 = cos(t617);
t607 = cos(t616);
t606 = cos(t615);
t605 = sin(t617);
t604 = sin(t616);
t603 = sin(t615);
t575 = -t584 * t626 + t587 * t629;
t574 = -t583 * t625 + t586 * t628;
t573 = -t582 * t624 + t585 * t627;
t569 = -t584 * t605 + t587 * t608;
t568 = -t583 * t604 + t586 * t607;
t567 = -t582 * t603 + t585 * t606;
t566 = t584 * t608 + t587 * t605;
t565 = t583 * t607 + t586 * t604;
t564 = t582 * t606 + t585 * t603;
t1 = [0, t631, t573 * t640 + t574 * t638 + t575 * t636, t631 * pkin(1), 0, t564 * t640 + t565 * t638 + t566 * t636 + (t564 * t644 + t565 * t643 + t566 * t642) * t630, t567 * t640 + t568 * t638 + t569 * t636 + (t567 * t644 + t568 * t643 + t569 * t642) * t630, -g(1); 0, t632, t573 * t641 + t574 * t639 + t575 * t637, t632 * pkin(1), 0, t564 * t641 + t565 * t639 + t566 * t637 + (t564 * t647 + t565 * t646 + t566 * t645) * t630, t567 * t641 + t568 * t639 + t569 * t637 + (t567 * t647 + t568 * t646 + t569 * t645) * t630, -g(2); 0, 0, 0, -3 * g(3), 0, 0, 0, -g(3);];
tau_reg  = t1;
