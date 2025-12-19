#pragma once

#include <QSlider>
#include <QVBoxLayout>
#include <QWidget>

#include <iostream>

/// A Qt widget for controlling a double value within a range [min, max].
/// The widget provides a similar interface as QSlider.
class DoubleSlider : public QWidget {
	Q_OBJECT

  private:
	void initialize() {
		m_intSlider->setRange(0, m_precision);
		auto* layout = new QVBoxLayout(this);
		layout->addWidget(m_intSlider);
		layout->setContentsMargins(0, 0, 0, 0);
		setSizePolicy(m_intSlider->sizePolicy());
		connect(m_intSlider, &QSlider::valueChanged, this, &DoubleSlider::handleIntValueChanged);
	}

  public:
	DoubleSlider(QWidget* parent = nullptr) : QWidget(parent) {
		m_intSlider = new QSlider(this);
		initialize();
	}

	DoubleSlider(Qt::Orientation orientation, QWidget* parent = nullptr) : QWidget(parent) {
		m_intSlider = new QSlider(orientation, this);
		initialize();
	}

	QSize sizeHint() const override {
		return m_intSlider->sizeHint();
	}

	void setOrientation(Qt::Orientation orientation) {
		m_intSlider->setOrientation(orientation);
		setSizePolicy(m_intSlider->sizePolicy());
	}

	void setMinimum(double min) {
		m_min = min;
	}

	void setMaximum(double max) {
		m_max = max;
	}

	void setRange(double min, double max) {
		m_min = min;
		m_max = max;
	}

	/// Set the number of discrete values (+1) that the slider uses.
	/// Example: a precision of 1 has 2 discrete steps: the minimum and maximum value.
	/// The default precision is 1000.
	void setPrecision(int precision) {
        double oldValue = value();
		m_precision = precision;
        m_intSlider->setMaximum(precision);
        setValue(oldValue);
	}

	void setValue(double val) {
		int sliderVal = static_cast<int>(((val - m_min) / (m_max - m_min)) * m_precision);
		m_intSlider->setValue(sliderVal);
	}

	[[nodiscard]] double value() const {
		return m_min + (m_intSlider->value() / static_cast<double>(m_precision)) * (m_max - m_min);
	}

    [[nodiscard]] double minimum() const {
        return m_min;
    }

    [[nodiscard]] double maximum() const {
        return m_max;
    }

  signals:
	void valueChanged(double value);

  private slots:
	void handleIntValueChanged(int) {
		emit valueChanged(value());
	}

  private:
	double m_min = 0.0;
	double m_max = 1.0;
	QSlider* m_intSlider;
	int m_precision = 1000;
};
