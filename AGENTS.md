# AGENTS.md — 전자공학 시각화 튜터

이 저장소에서 작업하는 에이전트(Hermes 등)가 읽고 따라야 할 규칙입니다.

## 이 저장소가 무엇인가

김준형의 전자공학 학습용 정적 웹 도구 모음입니다.
서버·빌드·외부 라이브러리 없이 **HTML 파일 하나가 곧 도구 하나**입니다.

| 파일 | 과목 | 내용 |
|---|---|---|
| `mosfet.html` | 반도체소자 | n-MOSFET 컷오프/선형/포화, 출력·전달 특성, 채널 에너지(핀치오프), g_m·r_o |
| `pmos-cmos.html` | 반도체소자 | p-MOSFET 부호규약, CMOS 인버터 VTC, 스위칭 문턱 V_M |
| `mos-cap.html` | 반도체소자 | MOS 커패시터 축적/평탄/공핍/반전, C-V, 게이트-산화막 개략 |
| `schottky-band.html` | 반도체소자 | Schottky 밴드 다이어그램 (평형/역바이어스) |
| `three-phase.html` | 회로이론 | Y-Y·Y-Δ·Δ-Y·Δ-Δ, 평형/불평형, 4선식/3선식, 페이저도, 풀이 과정 |
| `circuit-editor.html` | 회로이론 | R·L·C 회로를 직접 그리고 테브난/노턴(교류 복소) 구하기 |
| `thevenin.html` | 회로이론 | 값 편집형 테브난/노턴 단계 확인 |
| `ohms-law.html` | 회로이론 | 옴의 법칙 V=IR, 직렬·병렬·분압·전력, 예제 4종 단계별 풀이 |
| `maxwell.html` | 전자기학·전자장 | 맥스웰 방정식 4종(적분·미분형)과 각 식의 수치 예시 |
| `convolution.html` | 신호받시스템 | 컨볼루션 뒤집기·이동·곱하기·적분, 해석해 대조 |
| `eigen.html` | 선형대수 | 2×2 고유값·고유벡터, 단위원→타원, 복소 고유값 |
| `gauss.html` | 전자장 | 가우스 법칙 — 점·선·평면·구 대칭, 풀이 과정 |
| `pipeline.html` | 컴퓨터구조 | 5단 파이프라인 RAW 해저드, 스톨·포워딩, CPI |
| `index.html` | — | 도구 목록 페이지 (새 도구를 추가하면 여기에 카드를 넣는다) |

## 절대 규칙

1. **단일 파일 유지.** 각 도구는 `무언가.html` 하나로 완결된다. 외부 CDN, npm, 빌드 단계 금지.
2. **검증 없이 "된다"고 말하지 않는다.** 아래 검증 스크립트를 실제로 돌리고 종료 코드를 확인한다.
3. **숫자는 독립 계산과 대조한다.** 같은 코드로 낸 값을 같은 코드로 확인하는 것은 검증이 아니다.
4. **기존 페이지를 깨지 않는다.** `index.html`의 카드 링크가 모두 200으로 열려야 한다.
5. **개인정보·토큰을 커밋하지 않는다.** 이 저장소는 공개다.

## 검증 스크립트

```bash
python verify_mosfet.py        # n-MOSFET 4개 동작점, 경계 연속성, g_m/r_o 수치미분
python verify_pmos_cmos.py     # pMOS 3개 동작점, CMOS VTC 7점, V_M, 전류평형
python verify_mosc.py          # MOS-C 7항목
python verify_tyu91.py         # 교재 Schottky 예제 6개 값 + 불변식 5개
python verify_3phase.py        # 3상 4가지 결선 기준값
python verify_thevenin.py      # R_th를 세 경로(테스트전원/I_sc/단계축소)로 대조
python verify_ac_thevenin.py   # R·L·C 교류 테브난, 복소 정합
python verify_convolution.py   # 컨볼루션: 해석해 대조·교환법칙·δ 항등식·격자 수렴
python verify_eigen.py         # 특성다항식 잔차·‖Av−λv‖·Vieta·거듭제곱법·회전각
python verify_gauss.py         # 플럭스 직접 적분=Q/ε₀·껍질 상쇄·MC Q_enc·선·평면 대조
python verify_pipeline.py      # 손계산 8케이스 스케줄·구조 항등식·해저드 원장부
python verify_ohms_law.py      # MNA 노달 해석 vs 닫힌형·전력 삼경로·예제 4종 답
python verify_maxwell.py       # 맥스웰 4종 수치 대조·쌍극자 ∮B·dA=0·I_d=I_cond·빛의 속도
```

모두 표준 라이브러리만 사용하고, **종료 코드 0 = 통과**입니다.
도구를 수정하면 해당 검증 스크립트의 기대값도 함께 고치고, 왜 바뀌었는지 커밋 메시지에 남긴다.

## 화면을 고칠 때 지켜온 설계

- **"한 단계씩 보여주기"가 이 프로젝트의 핵심 요구다.** 결과만 보여주면 실패로 본다.
- 계산에 쓴 중간값을 화면에 그대로 노출한다(그래야 결과와 과정이 어긋나지 않는다).
- 물리적으로 불가능한 그림을 그리지 않는다. 축을 압축했으면 **압축했다고 표시**한다.
  (예: Schottky 역바이어스에서 밴드 굽힘이 Eg의 4.7배 → 세로축 압축 배지)
- 모델의 한계(1차 square-law, 종속전원 미지원 등)를 화면이나 README에 적는다.
- 한국어 UI. 단위(Ω, mH, µF, V)를 라벨에 명시한다.

## 배포

`main` 브랜치에 커밋하면 GitHub Pages가 자동으로 다시 빌드한다(약 1분).
사이트: https://junuary-k.github.io/electronics-tutor/

```bash
git add -A
git commit -m "..."
git push
```

## 확인 방법 (배포 후)

```bash
curl -s -o /dev/null -w "%{http_code}\n" https://junuary-k.github.io/electronics-tutor/
```

각 페이지가 200인지 확인한다. 404면 Pages 빌드가 아직 안 끝났거나 파일명이 틀린 것이다.
